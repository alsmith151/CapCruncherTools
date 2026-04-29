use itertools::Itertools;
use polars::prelude::*;
use pyo3_polars::PyDataFrame;
use rayon::prelude::*;
use std::collections::HashMap;

pub fn count(df: DataFrame) -> PyDataFrame {
    let parent_ids = df
        .column("parent_id")
        .expect("couldnt extract parent_id column")
        .cast(&DataType::UInt64)
        .expect("Failed to cast parent_id into u64");
    let restriction_fragments = df
        .column("restriction_fragment")
        .expect("couldnt extract restriction_fragment column")
        .cast(&DataType::Int64)
        .expect("Failed to cast restriction_fragment into i64");

    let parent_ids = parent_ids.u64().expect("parent_id should be u64");
    let restriction_fragments = restriction_fragments
        .i64()
        .expect("restriction_fragment should be i64");

    let mut fragments_by_parent: HashMap<u64, Vec<i64>> = HashMap::new();
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
                    *a.entry(k).or_insert(0) += v;
                }
                a
            },
        )
        .reduce(HashMap::new, |mut a, b| {
            for (k, v) in b {
                *a.entry(k).or_insert(0) += v;
            }
            a
        });

    let df_counts = DataFrame::new(
        interaction_counts.len(),
        vec![
            Column::new(
                "bin1_id".into(),
                interaction_counts
                    .keys()
                    .map(|(a, _)| *a)
                    .collect::<Vec<_>>(),
            ),
            Column::new(
                "bin2_id".into(),
                interaction_counts
                    .keys()
                    .map(|(_, b)| *b)
                    .collect::<Vec<_>>(),
            ),
            Column::new(
                "count".into(),
                interaction_counts.values().map(|v| *v).collect::<Vec<_>>(),
            ),
        ],
    )
    .expect("couldnt create dataframe");

    PyDataFrame(df_counts)
}
