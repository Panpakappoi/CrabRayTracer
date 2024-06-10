pub mod util{
    const EPSILON : f32 = 0.00001;
    pub fn equality(a:f32, b:f32) -> bool{
        let difference = a - b;
        if difference.abs() <= EPSILON{
            true
        }else{
            false
        }
    }
}