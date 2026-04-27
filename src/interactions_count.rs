use itertools::Itertools;
use polars::prelude::*;
use pyo3_polars::PyDataFrame;
use rayon::prelude::*;
use std::collections::HashMap;

pub fn count(df: DataFrame) -> PyDataFrame {
    let parent_ids = df
        .column("parent_id")
        .expect("couldnt extract parent_id column")
        .i64()
        .expect("parent_id should be i64");
    let restriction_fragments = df
        .column("restriction_fragment")
        .expect("couldnt extract restriction_fragment column")
        .i64()
        .expect("restriction_fragment should be i64");

    let mut fragments_by_parent: HashMap<i64, Vec<i64>> = HashMap::new();
    for (parent_id, restriction_fragment) in parent_ids
        .into_iter()
        .zip(restriction_fragments.into_iter())
    {
        if let (Some(parent_id), Some(restriction_fragment)) = (parent_id, restriction_fragment) {
            fragments_by_parent
                .entry(parent_id)
                .or_default()
                .push(restriction_fragment);
        }
    }

    let interaction_counts: HashMap<(i64, i64), i32> = fragments_by_parent
        .into_par_iter()
        .map(|(_, rf)| {
            let mut rf_combs = HashMap::new();
            for comb in rf.iter().combinations(2) {
                let (a, b) = vec![*comb[0], *comb[1]]
                    .into_iter()
                    .sorted()
                    .collect_tuple()
                    .unwrap();
                *rf_combs.entry((a, b)).or_insert(0) += 1;
            }

            rf_combs
        })
        .fold(
            HashMap::new,
            |mut a: HashMap<(i64, i64), i32>, b: HashMap<(i64, i64), i32>| {
                for (k, v) in b {
                    *a.entry(k).or_insert(1) += v;
                }
                a
            },
        )
        .reduce(HashMap::new, |mut a, b| {
            for (k, v) in b {
                *a.entry(k).or_insert(1) += v;
            }
            a
        });

    let height = interaction_counts.len();
    let df_counts = DataFrame::new(
        height,
        vec![
            Series::new(
                "bin1_id".into(),
                interaction_counts
                    .keys()
                    .map(|(a, _)| *a)
                    .collect::<Vec<_>>(),
            )
            .into(),
            Series::new(
                "bin2_id".into(),
                interaction_counts
                    .keys()
                    .map(|(_, b)| *b)
                    .collect::<Vec<_>>(),
            )
            .into(),
            Series::new(
                "count".into(),
                interaction_counts.values().map(|v| *v).collect::<Vec<_>>(),
            )
            .into(),
        ],
    )
    .expect("couldnt create dataframe");

    PyDataFrame(df_counts)
}
