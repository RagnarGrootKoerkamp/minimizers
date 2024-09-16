#![allow(dead_code)]

mod deinterleave;
mod gather;
mod table_lookup;

pub use deinterleave::deinterleave;
pub use gather::gather;
pub use table_lookup::table_lookup;
