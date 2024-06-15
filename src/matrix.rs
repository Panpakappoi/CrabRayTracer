use std::ops::{Mul};
use crate::tuple;
use crate::tuple::Tuple;
use crate::util::util::get_matrix;

pub fn translation(x: f32, y: f32, z: f32) -> Mat {
    Mat::new(MatDim::FourDim(
        vec![(1., 0.,0.,0.)],
        vec![(0.,1.,0.,0.)],
        vec![(0.,0.,1.,0.)],
        vec![(x,y,z,1.)]
    ))
}
pub fn inverse_translation(x: f32, y:f32, z:f32) -> Mat{
    Mat::new(MatDim::FourDim(
        vec![(1.,0.,0.,0.)],
        vec![(0.,1.,0.,0.)],
        vec![(0.,0.,1.,0.)],
        vec![(x * -1., y *-1., z *-1., 1.)]
    ))
}

pub fn scaling(x:f32, y: f32, z: f32) -> Mat{
    Mat::new(MatDim::FourDim(
        vec![(x, 0.,0.,0.)],
        vec![(0.,y,0.,0.)],
        vec![(0.,0.,z,0.)],
        vec![(0.,0.,0.,1.)]
    ))
}

pub fn shearing(x_y : f32, x_z : f32,
                y_x: f32, y_z : f32,
                z_x: f32, z_y : f32) -> Mat{
    Mat::new(MatDim::FourDim(
        vec![(1., y_x, z_x, 0.)],
        vec![(x_y, 1., z_y, 0.)],
        vec![(x_z, y_z, 1., 0.)],
        vec![(0.,0.,0.,1.)]
    ))
}
pub fn rotation_x(radians: f32) -> Mat{
    let _a = Mat::new(MatDim::FourDim(
        vec![(1.,0.,0.,0.)],
        vec![(0., radians.cos(), radians.sin(), 0.)],
        vec![(0. , (-1.* radians).sin(), radians.cos(), 0.)],
        vec![(0., 0., 0., 1.)]
    ));
    _a
}
pub fn rotation_y(radians: f32) -> Mat{
    let _a = Mat::new(MatDim::FourDim(
        vec![(radians.cos(),0.,(-1.0*radians).sin(),0.)],
        vec![(0., 1., 0., 0.)],
        vec![( radians.sin() , 0., radians.cos(), 0.)],
        vec![(0., 0., 0., 1.)]
    ));
    _a
}

pub fn rotation_z(radians:f32) -> Mat{
    let _a = Mat::new(MatDim::FourDim(
        vec![(radians.cos(), radians.sin(), 0., 0.)],
        vec![((-1. * radians).sin(), radians.cos(), 0., 0.)],
        vec![(0. , 0., 1., 0.)],
        vec![(0., 0., 0., 1.)]
    ));
    _a
}

#[derive(Debug, Clone)]
pub struct Mat {
    m_rows: usize,
    m_cols: usize,
    pub m_data : Vec<Vec<f32>>, // Im not happy about using this.
}

pub enum MatDim {
    // Im especially not happy with all the redundant memory allocations.
    // The problem is, is that tuple doesn't have an iterator, and that
    // the raw vec of floating points, makes it a little more difficult to index
    // when we're actually constructing the matrix.
    // This abomination is like, the compromise. But... As you can imagine,
    // Lots of computations, lots of allocations, probably not going to be a great time. // Revisit this i guess
    TwoDim(Vec<(f32, f32)>, Vec<(f32, f32)>),
    ThreeDim(Vec<(f32, f32, f32)>, Vec<(f32, f32, f32)>, Vec<(f32, f32, f32)>),
    FourDim(Vec<(f32, f32, f32, f32)>, Vec<(f32, f32, f32, f32)>,
            Vec<(f32,f32,f32,f32)>, Vec<(f32,f32,f32,f32)>)
}

impl Mul for Mat {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output {
        let mut matrix_data = vec![vec![0.;4], vec![0.;4], vec![0.;4], vec![0.;4]];
        for row in 0..self.m_rows{
            for col in 0..self.m_cols{
                matrix_data[row][col] = self.m_data[row][0] * rhs.m_data[0][col] +
                    self.m_data[row][1] * rhs.m_data[1][col] + self.m_data[row][2] *
                    rhs.m_data[2][col] + self.m_data[row][3] * rhs.m_data[3][col];
            }
        }
        Mat {m_rows: 4, m_cols:4, m_data : matrix_data}
    }
}
impl Mul<tuple::Tuple> for Mat {
    type Output = tuple::Tuple;
    // Reconsider here... The Tuple multiplication is fucked on rotation. Or its just..
    // With regard to the rotation matrix.. I wonder if we treat the radians as negative?
    fn mul(self, rhs: tuple::Tuple) -> Self::Output {
        let mut tuple_data:Vec<f32>= vec![] ;
        for a in self.m_data{
            let tuple_datum = a[0] * rhs.get_X() + a[1]*rhs.get_Y() + a[2] * rhs.get_Z() + a[3] *rhs.get_W();
            tuple_data.push(tuple_datum);
        }
        Tuple::new( tuple_data[0], tuple_data[1],tuple_data[2],  tuple_data[3])
    }
}

impl PartialEq<Vec<Vec<f32>>> for Mat {
    fn eq(&self, rhs: &Vec<Vec<f32>>) -> bool {
        if self.m_rows != rhs.len() || self.m_cols != rhs.len() {
            return false;
        }
        for i in 0..self.m_rows{
            for j in 0..self.m_cols{
                if self.m_data[j][i] != rhs[j][i]{
                    false;
                }
            }
        }
        true
    }
}

impl PartialEq for Mat {
    fn eq(&self, rhs: &Self) -> bool {
        if self.m_rows != rhs.m_rows || self.m_cols != rhs.m_cols {
            return false;
        }
        for i in 0..self.m_rows{
            for j in 0..self.m_cols{
                if self.m_data[j][i] != rhs.m_data[j][i]{
                    false;
                }
            }
        }
        true
    }
    fn ne(&self, rhs: &Self) -> bool{
        if self.m_rows != rhs.m_rows || self.m_cols != rhs.m_cols {
            return true;
        }
        for i in 0..self.m_rows{
            for j in 0..self.m_cols{
                if self.m_data[j][i] != rhs.m_data[j][i]{
                    return true;
                }
            }
        }
        false
    }
}

impl Mat {
    pub fn new(dim: MatDim) -> Self {
        match dim {
            MatDim::TwoDim(col1, col2) => {
                let mut data = vec![vec![], vec![]];
                for (v1, v2) in col1.into_iter().chain(col2) {
                    data[0].push(v1);
                    data[1].push(v2);
                }
                Mat { m_rows: 2, m_cols: 2, m_data: data }
            },
            MatDim::ThreeDim(col1, col2, col3) => {
                let mut data = vec![vec![], vec![], vec![]];
                for (v1, v2, v3) in col1.into_iter().chain(col2).chain(col3) {
                    data[0].push(v1);
                    data[1].push(v2);
                    data[2].push(v3);
                }
                Mat { m_rows: 3, m_cols: 3, m_data: data }
            }
            MatDim::FourDim(col1, col2, col3, col4) => {
                let mut data = vec![vec![], vec![], vec![], vec![]];
                for (v1, v2, v3, v4) in col1.into_iter().chain(col2).chain(col3).chain(col4) {
                    data[0].push(v1);
                    data[1].push(v2);
                    data[2].push(v3);
                    data[3].push(v4);
                }
                Mat { m_rows: 4, m_cols: 4, m_data: data }
            }
        }
    }
    pub fn transpose(&mut self) {
        let rows = self.m_data.len();
        let cols = self.m_data[0].len();
        for i in 0..rows {
            for j in (i + 1)..cols {
                let tmp = self.m_data[i][j];
                self.m_data[i][j] = self.m_data[j][i];
                self.m_data[j][i] = tmp;
            }
        }
    }
    pub fn get_determinant(&self) -> f32 {
        let result = match self.m_rows {
            // Hard indexing but when you expand to higher row and col numbers then... you're going to
            // have to come up with a more general solution?
            2 => (self.m_data[0][0] * self.m_data[1][1]) - (self.m_data[1][0] * self.m_data[0][1]),
            3 => {
                let _cof_a = self.get_cofactor(0,0);
                let _cof_b = self.get_cofactor(0,1);
                let _cof_c = self.get_cofactor(0, 2);
                (_cof_a * self.m_data[0][0]) + (_cof_b * self.m_data[0][1]) + (_cof_c * self.m_data[0][2])
            }
            4 => {
                let _cof_a = self.get_cofactor(0,0);
                let _cof_b = self.get_cofactor(0,1);
                let _cof_c = self.get_cofactor(0,2);
                let _cof_d = self.get_cofactor(0,3);
                (_cof_a * self.m_data[0][0]) +
                    (_cof_b * self.m_data[0][1]) +
                    (_cof_c * self.m_data[0][2]) +
                    (_cof_d * self.m_data[0][3])
            }
            _ => 0.0,
        };
        result
    }
    pub fn get_submatrix(&self, o_row: usize, o_col: usize) -> Vec<Vec<f32>> {
        self.m_data.iter()
            .enumerate()
            .filter_map(|(i, row)| {
                if i == o_row {
                    None
                } else {
                    Some(
                        row.iter()
                            .enumerate()
                            .filter_map(|(j, &val)| if j == o_col { None } else { Some(val) })
                            .collect::<Vec<f32>>()
                    )
                }
            })
            .collect()
    }
    pub fn get_minor(&self, o_row: usize, o_col: usize) -> f32{
        let _b = self.get_submatrix(o_row, o_col);
        let mut _b_sub = self.clone();
        _b_sub.m_cols = self.m_cols -1;
        _b_sub.m_rows = self.m_rows -1;
        _b_sub.m_data = _b;
        _b_sub.get_determinant()
    }
    pub fn get_cofactor(&self, o_row: usize, o_col: usize)->f32{
        let _minor = self.get_minor(o_row, o_col);
        if o_row + o_col % 2 == 1{
            return _minor * -1.
        } else {
            _minor
        }
    }
    pub fn is_invertible(&self) -> bool{
        if self.get_determinant() == 0.{
            false
        }else{
            true
        }
    }
    pub fn inverse(&self) -> Option<Self>{
        if !self.is_invertible(){
            None
        } else{
            let _det = self.get_determinant();
            let mut _b = self.clone();
            for i in 0..self.m_rows{
                for j in 0..self.m_cols{
                    let _c = self.get_cofactor(i,j);
                    _b.m_data[i][j] = _c / _det;
                }
            }
            Some(_b)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn construct_2d_mat() {
        let _2d_mat = Mat::new(MatDim::TwoDim(vec![(-3., 1.)], vec![(5., -2.)]));
        assert_eq!(_2d_mat.m_data[0][0], -3.);
        assert_eq!(_2d_mat.m_data[0][1], 5.);
        assert_eq!(_2d_mat.m_data[1][0], 1.);
        assert_eq!(_2d_mat.m_data[1][1], -2.);
    }

    #[test]
    fn construct_3d_mat() {
        let _3d_mat = Mat::new(MatDim::ThreeDim(vec![(-3., 1., 0.)], vec![(5., -2., 1.)], vec![(0., -7., 1.)]));
        assert_eq!(_3d_mat.m_data[0][0], -3.);
        assert_eq!(_3d_mat.m_data[1][1], -2.);
        assert_eq!(_3d_mat.m_data[2][2], 1.);
    }

    #[test]
    fn mat_equality_ident() {
        let _a =
            Mat::new(MatDim::FourDim(
                vec![(1., 5., 9., 5.)],
                vec![(2., 6., 8., 4.)],
                vec![(3., 7., 7., 3.)],
                vec![(4., 8., 6., 2.)]));
        let _b =
            Mat::new(MatDim::FourDim(vec![(1., 5., 9., 5.)],
                                     vec![(2., 6., 8., 4.)],
                                     vec![(3., 7., 7., 3.)],
                                     vec![(4., 8., 6., 2.)]));
        assert_eq!(_a == _b, true);
    }

    #[test]
    fn mat_equality_diff() {
        let _a =
            Mat::new(MatDim::FourDim(
                vec![(1., 5., 9., 5.)],
                vec![(2., 6., 8., 4.)],
                vec![(3., 7., 7., 3.)],
                vec![(4., 8., 6., 2.)]));
        let _b =
            Mat::new(MatDim::FourDim(
                vec![(2., 6., 8., 4.)],
                vec![(3., 7., 7., 3.)],
                vec![(4., 8., 6., 2.)],
                vec![(5., 9., 5., 1.)]));
        assert_eq!(_a != _b, true);
    }

    #[test]
    fn correct_mat_mult() {
        let _a =
            Mat::new(MatDim::FourDim(
                vec![(1., 5., 9., 5.)],
                vec![(2., 6., 8., 4.)],
                vec![(3., 7., 7., 3.)],
                vec![(4., 8., 6., 2.)]));
        let _b =
            Mat::new(MatDim::FourDim(vec![(-2., 3., 4., 1.)], vec![(1., 2., 3., 2.)], vec![(2., 1., 6., 7.)], vec![(3., -1., 5., 8.)]));
        assert_eq!(_a * _b, Mat::new(MatDim::FourDim(vec![(20., 44., 40., 16.)], vec![(22., 54., 58., 26.)], vec![(50., 114., 110., 46.)], vec![(48., 108., 102., 42.)])));
    }

    #[test]
    fn mat_tup_mult() {
        let _a =
            Mat::new(MatDim::FourDim(vec![(1., 2., 8., 0.)], vec![(2., 4., 6., 0.)], vec![(3., 4., 4., 0.)], vec![(4., 2., 1., 1.)]));
        let _b = Tuple::new(1., 2., 3., 1.);
        assert_eq!(_a * _b, Tuple::new(18., 24., 33., 1.));
    }

    #[test]
    fn mat_ident_mult() {
        let _a =
            Mat::new(MatDim::FourDim(vec![(0., 1., 2., 4.)], vec![(1., 2., 4., 8.)], vec![(2., 4., 8., 16.)], vec![(4., 8., 16., 32.)]));
        let _ident =
            Mat::new(MatDim::FourDim(vec![(1., 0., 0., 0.)], vec![(0., 1., 0., 0.)], vec![(0., 0., 1., 0.)], vec![(0., 0., 0., 1.)]));
        assert_eq!(_a * _ident, Mat::new(MatDim::FourDim(vec![(0., 1., 2., 4.)], vec![(1., 2., 4., 8.)], vec![(2., 4., 8., 16.)], vec![(4., 8., 16., 32.)])));
    }

    #[test]
    fn tuple_ident_mult() {
        let _a = Tuple::new(1., 2., 3., 4.);
        let _ident =
            Mat::new(MatDim::FourDim(vec![(1., 0., 0., 0.)], vec![(0., 1., 0., 0.)], vec![(0., 0., 1., 0.)], vec![(0., 0., 0., 1.)]));
        assert_eq!(_ident * _a, Tuple::new(1., 2., 3., 4.));
    }

    #[test]
    fn transpose_test() {
        let mut _a =
            Mat::new(MatDim::FourDim(vec![(0., 9., 1., 0.)], vec![(9., 8., 8., 0.)], vec![(3., 0., 5., 5.)], vec![(0., 8., 3., 8.)]));
        _a.transpose();
        assert_eq!(_a,
                   Mat::new(MatDim::FourDim(vec![(0., 9., 3., 0.)], vec![(9., 8., 0., 8.)], vec![(1., 8., 5., 3.)], vec![(0., 0., 5., 8.)])));
    }

    #[test]
    fn transpose_ident() {
        let mut _ident =
            Mat::new(MatDim::FourDim(vec![(1., 0., 0., 0.)], vec![(0., 1., 0., 0.)], vec![(0., 0., 1., 0.)], vec![(0., 0., 0., 1.)]));
        _ident.transpose();
        assert_eq!(_ident,
                   Mat::new(MatDim::FourDim(vec![(1., 0., 0., 0.)], vec![(0., 1., 0., 0.)], vec![(0., 0., 1., 0.)], vec![(0., 0., 0., 1.)])));
    }

    #[test]
    fn test_determinant() {
        let _a = Mat::new(MatDim::TwoDim(vec![(1., -3.)], vec![(5., 2.)]));
        let det = _a.get_determinant();
        assert_eq!(det, 17.)
    }

    #[test]
    fn get_sub_mat_test_3dim() {
        let _a = Mat::new(MatDim::ThreeDim(vec![(1., -3., 0.)], vec![(5., 2., 6.)], vec![(0., 7., -3.)]));
        let _a_sub = _a.get_submatrix(0, 2);
        println!("{:?}", _a_sub);
        assert_eq!(Mat::new(MatDim::TwoDim(vec![(-3., 0.)], vec![(2., 6.)])), _a_sub);
    }

    #[test]
    fn get_sub_mat_test_4dim() {
        let _a = Mat::new(MatDim::FourDim(vec![(-6., -8., -1., -7.)], vec![(1., 5., 0., 1.)], vec![(1., 8., 8., -1.)], vec![(6., 6., 2., 1.)]));
        let _a_sub = _a.get_submatrix(2, 1);
        let mut _b = _a.clone();
        _b.m_cols = 3;
        _b.m_rows = 3;
        _b.m_data = _a_sub;
        assert_eq!(Mat::new(MatDim::ThreeDim(vec![(-6., -8., -7.)], vec![(1., 8., -1.)], vec![(6., 6., 1.)])), _b);
    }

    #[test]
    fn minor_test() {
        let _a = Mat::new(MatDim::ThreeDim(vec![(3., 2., 6.)], vec![(5., -1., -1.)], vec![(0., -7., 5.)]));
        let _a_minor = _a.get_minor(1, 0);
        assert_eq!(_a_minor, 25.);
    }

    #[test]
    fn cofactor_test() {
        let _a = Mat::new(MatDim::ThreeDim(vec![(3., 2., 6.)], vec![(5., -1., -1.)], vec![(0., -7., 5.)]));
        let _a_cofactor = _a.get_cofactor(0, 0);
        let _a_cofactor_b = _a.get_cofactor(1, 0);
        assert_eq!(_a_cofactor, -12.);
        assert_eq!(_a_cofactor_b, -25.);
    }

    #[test]
    fn determinant_3d_test() {
        let _a =
            Mat::new(MatDim::ThreeDim(vec![(1., -5., 2.)], vec![(2., 8., 6.)], vec![(6., -4., 4.)]));
        let _a_det = _a.get_determinant();
        assert_eq!(_a_det, -196.);
    }

    #[test]
    fn determinant_4d_test() {
        let _a =
            Mat::new(MatDim::FourDim(
                vec![(-2., -3., 1., -6.)],
                vec![(-8., 1., 2., 7.)],
                vec![(3., 7., -9., 7.)],
                vec![(5., 3., 6., -9.)]));
        let _a_det = _a.get_determinant();
        assert_eq!(_a_det, -4071.);
    }

    #[test]
    fn invertible_test_succ() {
        let _a =
            Mat::new(MatDim::FourDim(
                vec![(6., 5., 4., 9.)],
                vec![(4., 5., -9., 1.)],
                vec![(4., 7., 3., 7.)],
                vec![(4., 6., -7., -6.)]));
        let b_invertible = _a.is_invertible();
        assert_eq!(b_invertible, true);
    }

    #[test]
    fn invertible_test_fail() {
        let _a =
            Mat::new(MatDim::FourDim(
                vec![(-4., 9., 0., 0.)],
                vec![(2., 6., -5., 0.)],
                vec![(-2., 2., 1., 0.)],
                vec![(-3., 6., -5., 0.)]));
        let b_invertible = _a.is_invertible();
        assert_eq!(b_invertible, false);
    }

    #[test]
    fn inverse_mat_test() {
        let _a =
            Mat::new(
                MatDim::FourDim(
                    vec![(8., 7., -6., -3.)],
                    vec![(-5., 5., 0., 0.)],
                    vec![(9., 6., 9., 9.)],
                    vec![(2., 1., 6., -4.)]));
        let _b = _a.inverse();
        let _c = get_matrix(_b);
        assert_eq!(_c,
                   Mat::new(MatDim::FourDim(vec![(-0.15385, -0.07692, 0.35897, -0.69231)],
                                            vec![(-0.15385, 0.1238, 0.35897, -0.69231)],
                                            vec![(-0.28205, 0.02564, 0.43590, -0.76923)],
                                            vec![(-0.53846, 0.03077, 0.92308, -1.92308)])));
    }

    #[test]
    fn inverse_mat_test_b() {
        let _a =
            Mat::new(MatDim::FourDim(
                vec![(9., -5., -4., -7.)],
                vec![(3., -2., 9., 6.)],
                vec![(0., -6., 6., 6.)],
                vec![(9., -3., 4., 2.)]));
        let _b = _a.inverse();
        let _c = get_matrix(_b);
        assert_eq!(_c,
                   Mat::new(MatDim::FourDim(vec![(-0.04074, -0.07778, -0.02901, 0.17778)],
                                            vec![(-0.07778, 0.03333, 0.36667, -0.33333)],
                                            vec![(-0.02901, -0.14630, -0.10926, 0.12963)],
                                            vec![(0.17778, 0.06667, -0.26667, 0.3333)])));
    }

    #[test]
    fn mat_product_inverse() {
        let _a = Mat::new(MatDim::FourDim(
            vec![(3., 3., -4., -6.)],
            vec![(-9., -8., 4., 5.)],
            vec![(7., 2., 4., -1.)],
            vec![(3., -9., 1., 1.)]
        ));
        let _b = Mat::new(MatDim::FourDim(
            vec![(8., 3., 7., 6.)],
            vec![(2., -1., 0., -2.)],
            vec![(2., 7., 5., 0.)],
            vec![(2., 0., 4., 5.)]
        ));
        // So... I mean ideally we implement traits for reference. Cause otherwise we're doing
        // alot of memory allocations. That's the thing with this entire project
        // Way too much allocations.
        let _c = _b.clone() * _a.clone();
        let _d = _c * get_matrix(_b.inverse());
        assert_eq!(_d, _a)
    }

    #[test]
    fn is_inv_trans_trans_inv() {
        let mut _a = Mat::new(
            MatDim::FourDim(
                vec![(9., -5., -4., -7.)],
                vec![(3., -2., 9., 6.)],
                vec![(0., -6., 6., 6.)],
                vec![(9., -3., 4., 2.)]));
        let mut _inv_b = _a.clone();
        _inv_b.transpose();
        let _inv_c = get_matrix(_inv_b.inverse());
        let mut _inv_a = get_matrix(_a.inverse());
        _inv_a.transpose();
        assert_eq!(_inv_a, _inv_c);
        // Proves that inverse transpose of matrix is transpose inverse of matrix
    }

    #[test]
    fn mult_translation() {
        let _transform = translation(5., -3., 2.);
        let _point = Tuple::new_point(-3., 4., 5.);
        let _transformed_point = _transform * _point;
        assert_eq!(_transformed_point, Tuple::new_point(2., 1., 7.));
    }

    #[test]
    fn mult_inverse_translation() {
        // There's something wrong with taking the translation matrix
        // and inverting it. It doesn't switch the signs of the x,y,z values
        // as is supposed to happen.
        // I think has something to do with cofactors. Which is strange
        // because the previous unit tests for matrix inversions worked fined.
        // I'm not sure.
        let mut _transform = inverse_translation(5., -3., 2.);
        println!("{:?}", _transform);
        let _point = Tuple::new_point(-3., 4., 5.);
        let _inv_point = _transform * _point;
        println!("{:?}", _inv_point);
        assert_eq!(_inv_point, Tuple::new_point(-8., 7., 3.));
    }

    #[test]
    fn mult_translation_vector() {
        let _transform = translation(5., 3., 2.);
        let _v = Tuple::new_vector(-3., 4., 5.);
        assert_eq!(_transform * _v, _v);
    }

    #[test]
    fn scale_point_test() {
        let _scaling = scaling(2., 3., 4.);
        let _p = Tuple::new_point(-4., 6., 8.);
        assert_eq!(_scaling * _p, Tuple::new_point(-8., 18., 32.));
    }

    #[test]
    fn scale_point_vector() {
        let _scaling = scaling(2., 3., 4.);
        let _p = Tuple::new_vector(-4., 6., 8.);
        assert_eq!(_scaling * _p, Tuple::new_vector(-8., 18., 32.));
    }

    #[test]
    fn scale_vector_inverse() {
        let _scaling = scaling(2., 3., 4.);
        let _scale_inv = get_matrix(_scaling.inverse());
        let _v = Tuple::new_vector(-4., 6., 8.);
        assert_eq!(_scale_inv * _v, Tuple::new_vector(-2., 2., 2.));
    }

    #[test]
    fn reflection_test() {
        let _scaling = scaling(-1., 1., 1.);
        let _p = Tuple::new_point(2., 3., 4.);
        assert_eq!(_scaling * _p, Tuple::new_point(-2., 3., 4.));
    }

    #[test]
    fn rotate_x(){
        let half_q = rotation_x(std::f32::consts::FRAC_PI_4);
        let full_q = rotation_x(std::f32::consts::FRAC_PI_2);
        let _p = Tuple::new_point(0.,1.,0.);
        assert_eq!(half_q * _p, Tuple::new_point(0., 2.0_f32.sqrt()/2., 2.0_f32.sqrt()/2.));
        assert_eq!(full_q * _p, Tuple::new_point(0.,0.,1.));
    }
    #[test]
    fn inverse_rotate_x(){
        let half_q = rotation_x(std::f32::consts::FRAC_PI_4).inverse();
        let _inv_half_q = get_matrix(half_q);
        let _p = Tuple::new_point(0.,1.,0.);
        assert_eq!(_inv_half_q * _p, Tuple::new_point(0., 2.0_f32.sqrt()/2., -2.0_f32.sqrt()/2.));
    }
    #[test]
    fn rotate_y(){
        let half_q = rotation_y(std::f32::consts::FRAC_PI_4);
        let full_q = rotation_y(std::f32::consts::FRAC_PI_2);
        let _p = Tuple::new_point(0.,0.,1.);
        assert_eq!(half_q * _p, Tuple::new_point( 2.0_f32.sqrt()/2.,0., 2.0_f32.sqrt()/2.));
        assert_eq!(full_q * _p, Tuple::new_point(1.,0.,0.));
    }
    #[test]
    fn rotate_z() {
        // So the whole sign system is still fucked. Uh...
        let half_q = rotation_z(std::f32::consts::FRAC_PI_4);
        let full_q = rotation_z(std::f32::consts::FRAC_PI_2);
        println!("{:?}", full_q);
        let _p = Tuple::new_point(0., 1., 0.);
        let _q = half_q.clone() * _p;
        println!("{:?}", _q);
        assert_eq!(half_q * _p, Tuple::new_point(-(2.0_f32.sqrt() / 2.), (2.0_f32.sqrt() / 2.), 0.));
        println!("{:?}", _p);
        assert_eq!(full_q * _p, Tuple::new_point(-1., 0., 0.));
    }
    #[test]
    fn shear_x_z_test(){
        let _trans = shearing(0.,1.,0.,0.,0.,0.);
        let _p = Tuple::new_point(2.,3.,4.);
        assert_eq!(_trans * _p, Tuple::new_point(6.,3.,4.));
    }
    #[test]
    fn shear_y_x_test(){
        let _trans = shearing(0.,0.,1.,0.,0.,0.);
        let _p = Tuple::new_point(2.,3.,4.);
        assert_eq!(_trans * _p, Tuple::new_point(2.,5.,4.));
    }
    #[test]
    fn shear_y_z_test(){
        let _trans = shearing(0., 0., 0., 1., 0., 0.);
        let _p = Tuple::new_point(2.,3.,4.);
        assert_eq!(_trans* _p, Tuple::new_point(2.,7.,4.));
    }
    #[test]
    fn shear_z_x_test(){
        let _trans = shearing(0.,0.,0.,0.,1.,0.);
        let _p = Tuple::new_point(2.,3.,4.);
        assert_eq!(_trans * _p , Tuple::new_point(2.,3.,6.));
    }
    #[test]
    fn shear_z_y_test(){
        let _trans = shearing(0.,0.,0.,0.,0.,1.);
        let _p = Tuple::new_point(2.,3.,4.);
        assert_eq!(_trans * _p, Tuple::new_point(2.,3.,7.));
    }
    #[test]
    fn mat_linear(){
        let _p = Tuple::new_point(1.,0.,1.);
        let _rot = rotation_x(std::f32::consts::FRAC_PI_2);
        let _scale = scaling(5.,5.,5.);
        let _trans = translation(10., 5.,7.);
        let _p2 = _rot* _p;
        assert_eq!(_p2, Tuple::new_point(1.,-1.,0.));
        let _p3 = _scale * _p2;
        assert_eq!(_p3, Tuple::new_point(5.,-5.,0.));
        let _p4 = _trans * _p3;
        assert_eq!(_p4, Tuple::new_point(15.,0.,7.));
    }
    #[test]
    fn mat_chain(){
        let _p = Tuple::new_point(1.,0.,1.);
        let _trans = translation(10.,5.,7.) *
            scaling(5.,5.,5.)*
            rotation_x(std::f32::consts::FRAC_PI_2) *
            _p;
        assert_eq!(_trans, Tuple::new_point(15.,0.,7.));
    }

}
