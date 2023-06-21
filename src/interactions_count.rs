use log::{debug, info, warn};
use polars::prelude::*;


struct ArgsPreProcessing{
    remove_exclusions: bool,
    remove_viewpoints: bool,
    subsample: bool,
    subsample_frac: Option<f64>,
}

impl ArgsPreProcessing{
}



fn remove_exclusions(df: DataFrame) -> Result<DataFrame> {
    df.filter(col("exclusion").ne("viewpoint"))
}

fn remove_viewpoints(df: DataFrame) -> Result<DataFrame>{
    df.filter(col("capture_count").eq(0))
}

fn subsample(df: DataFrame, frac: f64) -> Result<DataFrame>{
    df.sample_frac(frac)
}




fn get_viewpoint_from_parquet(path: &str, viewpoint: &str, args: ArgsPreProcessing) -> Result<()> {

    let df_lazy = LazyFrame::scan_parquet(path, Default::default())?
                     .filter(col("viewpoint").eq(viewpoint))
                     .select(&[col("viewpoint"), col("restriction_fragment")]);
    
    if args.remove_exclusions {
       let df_lazy = remove_exclusions(df_lazy)
    };






}
