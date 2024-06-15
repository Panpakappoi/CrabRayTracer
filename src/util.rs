use crate::matrix::Mat as Matrix;
pub mod util{
    use crate::matrix::Mat;
    use crate::matrix::MatDim;
    use crate::tuple::Tuple;
    const EPSILON : f32 = 0.00001;
    pub fn equality(a:f32, b:f32) -> bool{
        let difference = a - b;
        if difference.abs() <= EPSILON{
            true
        }else{
            false
        }
    }
    pub fn get_matrix(a: Option<Mat>) -> Mat {
        let b = a.unwrap_or(Mat::new(MatDim::FourDim(vec![(1., 0., 0., 0.)],
                                                     vec![(0.,1.,0.,0.)], vec![(0.,0.,1.,0.)], vec![(0.,0.,0.,1.)])));
        b
    }


}