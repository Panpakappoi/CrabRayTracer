
use std::ops::{Add, Sub, Neg, Mul, Div};

const EPSILON : f64 = 0.00001;
pub fn equality (a:f64, b:f64) -> bool { 
    let difference = a - b;
    if difference.abs() < EPSILON {
        true
    } else {
        false
    }
}


#[derive(Debug,PartialEq,Copy, Clone)]
pub struct Tuple {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub w: f64,
}

impl Tuple {
    pub fn new(x: f64, y: f64, z: f64, w:f64) -> Self
    {
        Tuple {
            x,
            y,
            z,
            w,
        }
    }

    pub fn new_point(x: f64, y: f64, z: f64) -> Self {
        Tuple::new(x, y, z, 1.0)
    }

    pub fn new_vector(x:f64, y: f64, z:f64) -> Self {
        Tuple::new(x,y,z, 0.0)
    }

    pub fn tuple_equality(&self, rhs: Tuple) -> bool {
        if equality(self.x, rhs.x) && equality(self.y, rhs.y) && equality(self.z, rhs.z) && equality(self.z, rhs.z)
        {
            true
        } else {
            false
        }
    }
    pub fn derive_magnitude(&self) -> f64 {
        let magnitude = (self.x * self.x + self.y*self.y + self. z * self.z + self. w * self.w).sqrt();
        magnitude
    }
    pub fn normalize(&self) -> Tuple {
        let magnitude = self.derive_magnitude();
        Tuple {
            x: self.x / magnitude,
            y: self.y / magnitude,
            z: self.z / magnitude,
            w: self.w / magnitude,
        }
    }
    pub fn dot_product(&self, rhs: Tuple) -> f64 {
        self.x * rhs.x + self.y * rhs.y + self.z * rhs.z + self.w * rhs.w
    }

    pub fn cross_product(&self, rhs: Tuple) -> Tuple {
        Tuple::new_vector(
            self.y * rhs.z - self.z * rhs.y,
            self.z * rhs.x - self.x * rhs.z,
            self.x *rhs.y - self.y * rhs.x
        ) 
    }
}

impl Add for Tuple {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        Tuple {x: self.x + rhs.x, y: self.y + rhs.y, z: self.z + rhs.z, w: self.w + rhs.w }
    } 
}

impl Sub for Tuple {
    type Output = Self;

    fn sub(self, rhs:Self) -> Self::Output {
        Tuple {x: self.x - rhs.x, y: self.y - rhs.y, z: self.z - rhs.z, w:self.w - rhs.w}
    } 
}
impl Neg for Tuple {
    type Output = Self; 
    fn neg(self) -> Self::Output {
        Tuple {x: -self.x, y: -self.y, z: -self.z, w: -self.w}
    }
}

impl Mul<f64> for Tuple {
    type Output = Tuple;
    fn mul(self, scalar: f64) -> Tuple {
        Tuple{x: self.x * scalar, y: self.y * scalar, z: self.z * scalar, w: self.w * scalar}
    }
}


trait DivBy<T> {
    type Output;
    fn div_by(self, scalar: T) -> Self::Output;
}

impl Div<f64> for Tuple {
    type Output = Tuple;
    fn div(self, scalar: f64) -> Tuple{
        Tuple{x:self.x / scalar, y: self.y / scalar, z: self.z / scalar, w: self.w / scalar }
    }
}

impl Div<i32> for Tuple {
    type Output = Tuple;
    fn div(self, scalar: i32) -> Tuple {
        Tuple{x:self.x / scalar as f64, y: self.y / scalar as f64, z: self.z / scalar as f64, w: self.w / scalar as f64}
    }
}

#[cfg(test)]
mod tests{
    use super::*;

    #[test]
    fn is_vector() {
        let _a = Tuple::new_point(4.3, -4.2, 3.1);
        assert_eq!(1.0, _a.w);
    }
    #[test] 
    fn is_point() {
        let _b = Tuple::new_vector(4.3, -4.2, 3.1);
        assert_eq!(0.0, _b.w);
    }
    #[test]
    fn addition_works(){
        let _a = Tuple::new_vector( 3.0, 2.0, 1.0);
        let _b = Tuple::new_point( -2.0, -1.0, 3.0);
        let result = Tuple {x: 1.0,y: 1.0,z: 4.0,w: 1.0};
        assert_eq!(_a + _b, result);
    }

    #[test]
    fn subtract_two_point(){
        let _a = Tuple::new_point(3.0, 2.0, 1.0);
        let _b = Tuple::new_point(5.0, 6.0, 7.0);
        let _result = Tuple::new_vector(-2.0, -4.0, -6.0);
        assert_eq!(_a - _b, _result);
    }

    #[test]
    fn subtract_vector_from_point() {
        let _p = Tuple::new_point(3.0, 2.0, 1.0);
        let _v = Tuple::new_vector(5.0, 6.0, 7.0);
        let _result = Tuple::new_point(-2.0, -4.0, -6.0);
        assert_eq!(_p-_v, _result);
    }

    #[test]
    fn subtraction_vector() {
        let _a = Tuple::new_vector(3.0, 2.0, 1.0);
        let _b = Tuple::new_vector(5.0, 6.0, 7.0);
        let result = Tuple{x:-2.0, y: -4.0, z:-6.0, w:0.0};
        assert_eq!(_a - _b, result);
    }
    #[test]
    fn negation_vector() {
        let _a = Tuple{x:5.0, y: 3.0, z: 6.0, w: 9.0};
        let negated_a = Tuple{x:-5.0, y:-3.0, z:-6.0, w:-9.0};
        assert_eq!(-_a, negated_a);
    }
    #[test]
    fn multiplying_a_tuple_by_scalar(){
        let _a = Tuple{x:1.0, y:-2.0, z:3.0,w:-4.0 };
        let scaled_a = Tuple{x:3.5, y:-7.0, z:10.5, w:-14.0};
        assert_eq!(_a * 3.5, scaled_a);
    }
    #[test]
    fn multiplying_a_tuple_by_fraction() {
        let _a = Tuple{x:1.0, y:-2.0, z:3.0, w:-4.0};
        let divided_a = Tuple{x:0.5, y:-1.0,z: 1.5, w:-2.0};
        assert_eq!(_a * 0.5, divided_a);
    }

    #[test]
    fn divide_by_scalar() {
        let _a = Tuple{x:1.0, y:-2.0, z:3.0, w:-4.0};
        let divided_a = Tuple{x:0.5, y:-1.0, z: 1.5, w:-2.0};
        assert_eq!(_a / 2, divided_a);
    }
    #[test]
    fn test_magnitude_a() {
        let _a = Tuple::new_vector(-1.0, -2.0, -3.0);
        let magnitude = 14_f64.sqrt();
        assert_eq!(_a.derive_magnitude(), magnitude);
    }
    #[test]
    fn test_magnitude_b() {
        let _b = Tuple::new_vector(0.0, 0.0, 1.0);
        let magnitude = 1_f64.sqrt();
        assert_eq!(_b.derive_magnitude(), magnitude);
    }
    #[test]
    fn test_normalize() {
        let _c = Tuple::new_vector(4.0, 0.0, 0.0);
        let _c_normalized = Tuple::new_vector(1.0, 0.0, 0.0);
        assert_eq!(_c.normalize(), _c_normalized);
    }
    #[test]
    fn test_normalize_b() {
        let _c = Tuple::new_vector(1.0,2.0,3.0);
        let _c_normalized = Tuple::new_vector(1.0/14.0_f64.sqrt(), 2.0/14.0_f64.sqrt(), 3.0/14.0_f64.sqrt());
        assert_eq!(_c.normalize(), _c_normalized);
    }
    #[test]
    fn test_noramlized_magnitude() {
        let mut _c = Tuple::new_vector(1.0, 2.0, 3.0);
        _c = _c.normalize();
        let magnitude = 1.0_f64;
        assert_eq!(_c.derive_magnitude(), magnitude);

    }

    #[test]
    fn test_dot_product() {
        let _c = Tuple::new_vector(1.0,2.0,3.0);
        let _d = Tuple::new_vector(2.0, 3.0, 4.0);
        let dp = 20.0_f64;
        assert_eq!(_c.dot_product(_d),dp);
    }
    #[test]
    fn test_cross_product() {
        let _c = Tuple::new_vector(1.0, 2.0, 3.0);
        let _d = Tuple::new_vector(2.0, 3.0, 4.0);
        let cp_1 = Tuple::new_vector(-1.0, 2.0, -1.0);
        assert_eq!(_c.cross_product(_d), cp_1, " Something went horribly wrong");
    }
    #[test]
    fn test_cross_product_b() {
        let _c = Tuple::new_vector(1.0, 2.0, 3.0);
        let _d = Tuple::new_vector(2.0, 3.0, 4.0);
        let cp_2 = Tuple::new_vector(1.0, -2.0, 1.0);
        assert_eq!(_d.cross_product(_c), cp_2, "Cross Product failed between {:?}, {:?}", _c, _d);
    }
}
