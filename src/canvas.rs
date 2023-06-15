use crate::util::util::equality;

use std::ops::{Add, Sub, Neg, Mul, Div};
#[derive(Debug,PartialEq,Copy, Clone)]
struct Color{
    red : f64,
    green : f64,
    blue : f64,
    
}

impl Color {
    fn new(red:f64, green:f64, blue:f64) -> Self{
        Color {
            red, 
            green,
            blue, 
            
        }
    }
    fn hadamard_product(&self, rhs: Color) -> Color {
        Color::new(self.red * rhs.red, self.green * rhs.green, self.blue * rhs.blue)
    }
}

impl Add for Color{
    type Output = Self;
    fn add(self, rhs:Self) -> Self::Output {
        Color::new(self.red + rhs.red, self.green + rhs.green, self.blue + rhs.blue)
    }
}

impl Sub for Color {
    type Output = Self;
    fn sub(self, rhs:Self) -> Self::Output {
        Color::new(self.red - rhs.red, self.green - rhs.green, self.blue - rhs.blue)
    }
}

impl Mul<f64> for Color {
    type Output = Self;
    fn mul(self, scalar: f64) -> Self::Output {
        Color::new(self.red * scalar, self.green * scalar, self.blue * scalar)
    }
}
impl Mul<Color> for Color{
    type Output = Self;
    fn mul(self, rhs: Self) -> Self::Output{
        Color::new(self.red * rhs.red, self.green * rhs.green, self.blue * rhs.blue)
    }
}

struct Canvas
{
    width: u32,
    height: u32,
    pixels: Vec<Color>,    
}

impl Canvas {
    fn new(width: u32, height: u32) -> Self {
        let pixels = vec![Color::new(0.0, 0.0, 0.0); (width * height) as usize];
        Canvas {
            width,
            height,
            pixels,
        }
    }

}


#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn Color_are_rgb_tuples() {
        let c = Color::new (-0.5, 0.4, 1.7);
        assert_eq!(c.red, -0.5);
        assert_eq!(c.green, 0.4);
        assert_eq!(c.blue, 1.7);
    }
    #[test]
    fn add_color_components() {
        let c1 = Color::new(0.9, 0.6, 0.75);
        let c2 = Color::new(0.7, 0.1, 0.25);
        let expected_result = Color::new (1.6, 0.7, 1.0);
        assert_eq!(c1+c2, expected_result);
    }
    #[test]
    fn subtract_color_components(){
        let c1 = Color::new(0.9, 0.6, 0.75);
        let c2 = Color::new(0.7, 0.1, 0.25);
        let actual_result = c1 - c2;
        let expected_result = Color::new (0.2, 0.5, 0.5);
        assert!(equality(actual_result.red, expected_result.red));
        assert!(equality(actual_result.green, expected_result.green));
        assert!(equality(actual_result.blue, expected_result.blue));
    }
    #[test]
    fn multiply_color_scalar(){
        let c = Color::new(0.2, 0.3, 0.4);
        let expected = Color::new(0.4, 0.6, 0.8);
        assert_eq!(c * 2.0, expected);
   } 
   #[test]
   fn multiply_color_by_color() {
    let c_1 = Color::new(1.0, 0.2, 0.4);
    let c_2 = Color::new(0.9, 1.0, 0.1);
    let actual = c_1 * c_2;
    let expected = Color::new(0.9, 0.2, 0.04) ;
    assert!(equality(actual.red, expected.red));
    assert!(equality(actual.green, expected.green));
    assert!(equality(actual.blue, expected.blue));
   }
   #[test]
   fn create_a_canvas() {
    let c = Canvas::new(10, 20);
    assert_eq!(c.width, 10);
    assert_eq!(c.height, 20);
    //assert each pixel is color (0,0,0)
    for pixel in c.pixels {
        assert_eq!(pixel, Color::new(0.0, 0.0, 0.0));
    }
   }
}
// impl Add for Color {
//     type Output = Self;

//     fn add(self, rhs: Self) -> Self::Output {
//         Tuple {x: self.x + rhs.x, y: self.y + rhs.y, z: self.z + rhs.z, w: self.w + rhs.w }
//     } 
