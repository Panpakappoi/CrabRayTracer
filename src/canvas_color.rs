

use std::fmt;
use crate::util::util::equality;
use std::ops::{Add, Sub, Mul, Div};
#[derive(Clone, Copy, Debug)]
pub struct Color{
    pub red : f32,
    pub green : f32,
    pub blue : f32
}

impl PartialEq for Color{
    fn eq(&self, other : &Self) -> bool{
        equality(self.red, other.red) && 
        equality(self.green, other.green) && 
        equality(self.blue, other.blue)
    }
}

impl Add for Color{
    type Output = Color;
    fn add(self, rhs: Color) -> Self{
        Color::new(self.red + rhs.red, self.green + rhs.green, self.blue + rhs.blue)
    }
}

impl Sub for Color{
    type Output = Color;
    fn sub(self, rhs: Color) -> Self{
        Color::new(self.red - rhs.red, self.green - rhs.green, self.blue - rhs.blue)
    }
}

impl Mul<f32> for Color{
    type Output = Color;
    fn mul(self, scalar:f32) -> Self{
        Color::new(self.red * scalar, self.green * scalar, self.blue * scalar)
    }
}

impl Mul<Color> for Color{
    type Output = Color;
    fn mul(self, rhs : Color) -> Self{
        Color::new(self.red * rhs.red, self.green* rhs.green, self.blue * rhs.blue)
    }
}

impl Color{
    pub fn new(red: f32, green : f32, blue : f32) -> Self{
        Color{
            red,
            green,
            blue
        }
    }
    pub fn hadamard_product(&self, rhs : Color)-> Self{
        Color::new(self.red * rhs.red, self.green * rhs.green, self.blue * rhs.blue)
    }

}
#[cfg(test)]
mod tests{ 
    use super::*;
    #[test]
    fn add_color_test(){
        let _c1 = Color::new(0.9, 0.6, 0.75);
        let _c2 = Color::new(0.7, 0.1, 0.25);
        assert_eq!(_c1 + _c2, Color::new(1.6,0.7,1.0));
    }
    #[test]
    fn sub_color_test(){
        let _c1 = Color::new(0.9, 0.6, 0.75);
        let _c2 = Color::new(0.7, 0.1, 0.25);
        assert_eq!(_c1 - _c2, Color::new(0.2, 0.5, 0.5));
    }
    #[test]
    fn mul_color_scalar_test(){
        let _c1 = Color::new(0.2, 0.3, 0.4);
        let scalar = 2.0;
        assert_eq!(_c1 * scalar, Color::new(0.4, 0.6, 0.8));
    }
}