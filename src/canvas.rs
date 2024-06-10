use crate::canvas_color::Color;
use std::fs::File;
use std::io::{Result, Write};
pub struct Canvas{
    width : u32,
    height : u32,
    pixels : Vec<Vec<Color>>
}

fn float_to_u8(norm_val:f32) -> u8{
    let mut n = norm_val;
    if n > 1.0{
        n = 1.0;
    } else if n < 0.0 {
        n = 0.0;
    }
    let scaled = n * 255.0;
    scaled.round() as u8
}

impl Canvas{
    pub fn new(height : u32, width : u32) -> Self{
        let pixels = vec![vec![Color::new(0.0,0.0,0.0); width as usize]; height as usize];
        Canvas{
            width,
            height,
            pixels
        }
    } // There really isn't a need for color to be f32?
    pub fn write_pixel(&mut self, c1: Color, w : usize, h: usize){
        self.pixels[h][w] = c1;
    }
    pub fn fill_canvas(&mut self, c1: Color){
        for i in 0..self.height as usize{
            for j in 0..self.width as usize{
                self.pixels[i][j] = c1;
            }
        }
    }
    pub fn write_ppm(self, filename: &str, width: usize, height : usize) -> Result<()>{
        let mut file = File::create(filename)?;
        writeln!(file, "P3")?;
        writeln!(file, "{} {}", width, height)?;
        writeln!(file, "255")?;
        let mut line_buf = String::new();
        for y in 0..height{
            for x in 0..width{
                let pixel_data = format!("{} {} {}", 
                float_to_u8(self.pixels[y][x].red) , float_to_u8(self.pixels[y][x].green),
                float_to_u8(self.pixels[y][x].blue));
                
                if line_buf.len() + pixel_data.len() + 1 > 70 {
                    writeln!(file, "{}", line_buf)?;
                    line_buf.clear();
                }

                if !line_buf.is_empty(){
                    line_buf.push(' ')
                }
                line_buf.push_str(&pixel_data);
            }
        }
        if !line_buf.is_empty(){
            writeln!(file, "{}", line_buf)?;
        }
        Ok(())
    }
    pub fn get_height(&self) -> u32{
        self.height
    }
    pub fn get_width(&self) -> u32{
        self.width
    }
}

#[cfg(test)]
mod tests{
    use super::*;
    #[test]
    fn create_canvas(){
        let canvas = Canvas::new(10,20);
        assert_eq!(canvas.height, 10);
        assert_eq!(canvas.width, 20);
        for row in canvas.pixels{
            for pixel in row{
                assert_eq!(pixel, Color::new(0.,0.,0.));
            }
        }
    }
    #[test]
    fn write_pixel_test(){
        let mut canvas = Canvas::new(10,20);
        canvas.write_pixel(Color::new(1.,0.,0.), 2,3);
        assert_eq!(canvas.pixels[3][2], Color::new(1.,0.,0.));
    }
    #[test]
    fn write_file_test(){
        let _canvas = Canvas::new(5,3);
    }
}

