use std::{fmt};

use crate::util::util::equality;
use std::ops::{Add, Sub, Neg, Mul, Div};


#[derive (Clone, Copy)]
pub struct Tuple{
    pub x: f32,
    pub y: f32,
    pub z: f32,
    pub w: f32
}

impl fmt::Debug for Tuple{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result{
        write!(f, "Tuple {{x: {}, y:{}, z:{}, w:{}}}", self.x, self.y, self.z, self.w)
    }
}

impl PartialEq for Tuple{
    fn eq(&self, other : &Self) -> bool{
        equality(self.x, other.x) && 
        equality(self.y, other.y) && 
        equality(self.z, other.z) && 
        equality(self.w, other.w)
    }
}

impl Add for Tuple{
    type Output = Tuple;

    fn add(self, rhs: Tuple) -> Self::Output {
        let mut _a = 
            Tuple::new(self.x + rhs.x, self.y + rhs.y, self.z + rhs.z, 0.0);
        if self.isPoint() || rhs.isPoint(){
            _a.set_W(1.0);
            _a
        }else {
            _a
        }
    }
}
impl Sub for Tuple{
    type Output = Tuple;
    fn sub(self, rhs: Self) -> Self::Output { 
            Tuple::new(self.x - rhs.x, self.y - rhs.y, self.z - rhs.z, self.w - rhs.w)
    }
}

impl Neg for Tuple{
    type Output = Tuple;
    fn neg(self) -> Self::Output{
        Tuple::new(self.x * -1.0, self.y * -1.0, self.z * -1.0, self.w * -1.0)
    } 
}

impl Mul<f32> for Tuple{
    type Output = Tuple;
    fn mul(self, scalar : f32) -> Self::Output{
        Tuple::new(self.x * scalar, self.y *scalar, self.z * scalar, self.w * scalar)
    }
}

impl Div<f32> for Tuple{
    type Output = Tuple;
    fn div(self, scalar: f32) -> Self::Output {
        Tuple::new(self.x / scalar, self.y /scalar, self.z / scalar, self.w / scalar )
    }
}

impl Tuple{
    // Constructor
    pub fn new(x :f32, y : f32, z: f32, w: f32) -> Tuple{
        Tuple{
            x,
            y,
            z,
            w,
        }
    }
    pub fn new_point(x:f32, y:f32, z:f32) -> Tuple{
        Tuple::new(x,y,z,1.0)
    }

    // Getters and Setters
    pub fn get_X(&self) -> f32{
        self.x
    }
    pub fn set_X(&mut self, x:f32){
        self.x = x;
    }
    pub fn get_Y(&self) -> f32{
        self.y
    }
    pub fn set_Y(&mut self, y: f32){
        self.y = y;
    }
    pub fn get_Z(&self) -> f32{
        self.z
    }
    pub fn set_Z(&mut self, z:f32){
        self. z = z;
    }
    pub fn get_W(&self) -> f32{
        self.w
    }
    pub fn set_W(&mut self, w: f32){
        self.w = w;
    }

    pub fn new_vector(x:f32, y:f32, z:f32) -> Tuple{
        Tuple::new(x,y,z,0.0)
    }
    pub fn isPoint(&self) -> bool {
        match self.get_W(){
            0.0 => false,
            1.0 => true,
            _ => false
        }
    } 
    pub fn isVector(&self) -> bool{
        match self.get_W(){
            0.0 => true,
            1.0 => false,
            _ => false
        }
    }
    pub fn isEqual(&self, rhs: Tuple) -> bool{
        if equality(self.x, rhs.x) && equality(self.y, rhs.y) && equality(self.z ,rhs.z) && equality(self.w, rhs.w){
            true
        } else{
            false
        }
    }

    pub fn magnitude(&self) -> f32{ 
        (self.x * self.x + self.y*self.y + self.z*self.z + self.w * self.w).sqrt() 
    }

    pub fn normalize(&self) -> Tuple{
        Tuple::new( self.x /self.magnitude(), self.y / self.magnitude(), self.z / self.magnitude(), self.w / self.magnitude())
    }

    pub fn dot_product(&self, rhs: Tuple) -> f32{
        self.x * rhs.x +
        self.y * rhs.y +
        self.z * rhs.z +
        self. w * rhs.w
    }
    pub fn cross_product(&self, rhs: &Tuple) -> Tuple{
        Tuple::new_vector(self.y * rhs.z - self.z * rhs.y,
                        self.z * rhs.x - self.x * rhs.z,
                        self.x * rhs.y - self.y * rhs.x)
    }
    

}

#[cfg(test)]
mod tests{
    use super::*;
    #[test]
    fn is_point_test(){
        let _a = Tuple::new_point(4.3, -4.2, 3.1);
        let isPoint = _a.isPoint();
        let isVector = _a.isVector();
        assert_eq!(isPoint, true);
        assert_eq!(isVector, false);
    }
    #[test]
    fn is_vector_test(){
        let _a  = Tuple::new_vector(4.3, -4.2, 3.1);
        let isPoint = _a.isPoint();
        let isVector = _a.isVector();
        assert_eq!(isPoint, false);
        assert_eq!(isVector, true);
    }
    #[test]
    fn add_two_tuple_test(){
        let _a1 = Tuple::new(3.0, -2.0, 5.0, 1.0);
        let _a2 = Tuple::new(-2.0, 3.0, 1.0, 0.0);
        let _b = Tuple::new(1.0, 1.0, 6.0, 1.0);
        let _c = _a1 + _a2;
        assert_eq!(_c, _b )
    }
    #[test]
    fn sub_two_point_test(){
        let _p1 = Tuple::new_point(3.0,2.0,1.0);
        let _p2 = Tuple::new_point(5.0,6.0,7.0);
        assert_eq!(_p1 - _p2, Tuple::new(-2.0, -4.0, -6.0, 0.0));
    }
    #[test]
    fn sub_point_vec_test(){
        let _p1 = Tuple::new_point(3.0, 2.0, 1.0);
        let _v1 = Tuple::new_vector(69.0, 42.0, 85.0);
        assert_eq!(_p1 - _v1, Tuple::new_point(-66.0, -40.0, -84.0));
    }
    #[test]
    fn sub_two_vec_test(){
        let _v1 = Tuple::new_vector(3.0, 2.0, 1.0);
        let _v2 = Tuple::new_vector(5.0, 6.0, 7.0);
        assert_eq!(_v1 - _v2 ,Tuple::new_vector(-2.0, -4.0, -6.0));
    }
    #[test]
    fn negation_test(){
        let _a1 = Tuple::new(1.0,2.0,3.0,4.0);
        assert_eq!(-_a1, Tuple::new(-1.0, -2.0, -3.0, -4.0));
    }

    #[test]
    fn scalar_mul_test(){
        let mut _a1 = Tuple::new(1.0, -2.0, 3.0, -4.0);
        _a1 = _a1 * 3.5;
        assert_eq!(_a1, Tuple::new(3.5, -7.0, 10.5, -14.0))
    }
    #[test]
    fn fractionated_mul_test(){
        let mut _a1 = Tuple::new(1.0, -2.0, 3.0, -4.0);
        _a1 = _a1 * 0.5;
        assert_eq!(_a1, Tuple::new(0.5, -1.0, 1.5, -2.0));
    }
    #[test]
    fn scalar_div_test(){
        let mut _a1 = Tuple::new(1.0,-2.0, 3.0, -4.0);
        _a1 = _a1 / 2.0;
        assert_eq!(_a1, Tuple::new(0.5, -1.0, 1.5, -2.0)); 
    }
    #[test]
    fn magnitude_1_test(){
        let _v = Tuple::new_vector(1.0, 0.0, 0.0);
        let mag = _v.magnitude();
        assert_eq!(mag, 1.0);
    }
    #[test]
    fn magnitude_2_test(){
        let _v = Tuple::new_vector(0.0, 1.0, 0.0);
        let mag = _v.magnitude();
        assert_eq!(mag,1.0);
    }
    #[test]
    fn magnitude_3_test(){
        let _v = Tuple::new_vector(0.0, 0.0, 1.0);
        let mag = _v.magnitude();
        assert_eq!(mag,1.0);
    }
    #[test]
    fn magnitude_4_test(){
        let _v = Tuple::new_vector(1.0, 2.0, 3.0);
        let mag = _v.magnitude();
        assert_eq!(mag,f32::sqrt(14.0));
    }

    #[test]
    fn magnitude_5_test(){
        let _v = Tuple::new_vector(-1.0, -2.0, -3.0);
        let mag = _v.magnitude();
        assert_eq!(mag,f32::sqrt(14.0));
    }
    #[test]
    fn normalize_1_test(){
        let _v = Tuple::new_vector(1.0,2.0,3.0);
        assert_eq!(_v.normalize(), Tuple::new_vector(0.26726, 0.53452, 0.80178));
    }
    #[test]
    fn normalize_2_test(){
       let _v = Tuple::new_vector(4.0,0.0,0.0);
        assert_eq!(_v.normalize(), Tuple::new_vector(1.0,0.0,0.0));
    }
    #[test]
    fn magnitude_of_norm_test(){
        let _v = Tuple::new_vector(1.0,2.0,3.0);
        let _v_norm = _v.normalize();
        assert_eq!(equality(_v_norm.magnitude(), 1.0_f32), true);
    }

    #[test]
    fn cross_two_vec_test(){
        let _a1 = Tuple::new_vector(1.0, 2.0, 3.0);
        let _a2 = Tuple::new_vector(2.0,3.0,4.0);
        let _a1_cross = _a1.cross_product(&_a2);
        let _a2_cross = _a2.cross_product(&_a1);
        assert_eq!(_a1_cross, Tuple::new_vector(-1.0, 2.0, -1.0));
        assert_eq!(_a2_cross, Tuple::new_vector(1.0, -2.0, 1.0));
    }




}
