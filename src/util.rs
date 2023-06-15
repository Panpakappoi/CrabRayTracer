pub mod util{
const EPSILON : f64 = 0.00001;
pub fn equality (a:f64, b:f64) -> bool { 
    let difference = a - b;
    if difference.abs() < EPSILON {
        true
    } else {
        false
    }
}}