use std::{error::Error, path::{PathBuf, Path}};

use anyhow::Ok;
use itertools::Itertools;
use log::{debug, info, warn};
use polars::prelude::*;
use pyo3::prelude::*;
use pyo3_polars::PyDataFrame;
use rayon::prelude::*;
use std::collections::HashMap;

pub fn count(df: DataFrame) -> PyDataFrame {
    // Group the dataframe by parent_id


    let interaction_counts: HashMap<(i64, i64), i32> = df
        .partition_by(vec!["parent_id"], false)
        .expect("couldnt partition by parent_id")
        .into_iter()
        .map(|df| {
            let ser = &df
                .select_series(["restriction_fragment"])
                .expect("couldnt extract restriction_fragment column")[0];
            let rf = ser.i64().expect("should be i64");

            let mut rf_combs = HashMap::new();
            for comb in rf.into_iter().combinations(2) {
                let (a, b) = (comb[0], comb[1]);

                match (a, b) {
                    (Some(a), Some(b)) => {
                        let (a, b) = vec![a, b].into_iter().sorted().collect_tuple().unwrap();
                        *rf_combs.entry((a, b)).or_insert(0) += 1;
                    }
                    _ => {}
                }
            }

            rf_combs
        })
        .fold(
            HashMap::new(),
            |mut a: HashMap<(i64, i64), i32>, b: HashMap<(i64, i64), i32>| {
                for (k, v) in b {
                    *a.entry(k).or_insert(1) += v;
                }
                a
            },
        );

    let df_counts = DataFrame::new(vec![
        Series::new(
            "bin1_id",
            interaction_counts
                .keys()
                .map(|(a, _)| *a)
                .collect::<Vec<_>>(),
        ),
        Series::new(
            "bin2_id",
            interaction_counts
                .keys()
                .map(|(_, b)| *b)
                .collect::<Vec<_>>(),
        ),
        Series::new(
            "count",
            interaction_counts.values().map(|v| *v).collect::<Vec<_>>(),
        ),
    ])
    .expect("couldnt create dataframe");

    PyDataFrame(df_counts)
}
