mod util;
mod projectile;
mod tuple;
mod canvas;
mod ray;
mod sphere;
mod canvas_color;
mod matrix;
use crate::projectile::Environment;
use crate::projectile::Projectile;
use crate::projectile::tick;
use crate::canvas_color::Color;
use crate::canvas::Canvas;

fn trace_proj(proj: Projectile, env : Environment, canvas : &mut Canvas, _c1:Color){
    let mut a = &proj;
    while(a.position.y > 0.){
       // a = tick(env, a);
        let height = if (a.position.y.round() as u32) < canvas.get_height() {
            a.position.y.round() as usize
        } else if a.position.y.round() as u32 > canvas.get_height() {
            canvas.get_height() as usize
        } else{
            0
        };
        let width = if (a.position.x.round() as u32) < canvas.get_width(){
            a.position.x.round() as usize
        } else if a.position.x.round() as u32 > canvas.get_width(){
            canvas.get_width() as usize
        } else{
            0
        };
        canvas.write_pixel(_c1, width, height);
    }

}

fn main() {

    let mut proj_canvas = Canvas::new(600, 1000);
    let mut a = Projectile::new(0.0,1.0,0.0, 1.1,1.8,0.0);
    a.velocity = a.velocity.normalize() * 11.25;
    let b = Environment::new(0., -0.1, 0., -0.01, 0., 0.);
    let mut counter = 0;
    let _c1 = Color::new(0.333, 1.0, 0.843);
    while(a.position.y > 0.){
        a = tick(&b, &a);
        let height = if (a.position.y.round() as u32) < proj_canvas.get_height() 
            && a.position.y.round() > 0.{
            a.position.y.round() as usize
        } else if a.position.y.round() as u32 > proj_canvas.get_height() {
            500
        } else{
            0
        };
        let width = if (a.position.x.round() as u32) < proj_canvas.get_width()
         && a.position.x.round() > 0.{
            a.position.x.round() as usize
        } else if a.position.x.round() as u32 > proj_canvas.get_width(){
            900
        } else{
            0
        };
        proj_canvas.write_pixel(_c1, width,
         height);
    }
    let succ_2 = proj_canvas.write_ppm("cannon.ppm", 1000, 550);
    

    let mut canva_1 = Canvas::new(2,10);
    let _c4 = Color::new(1.0,0.8,0.6);
    canva_1.fill_canvas(_c4);
    let _succ1 = canva_1.write_ppm("fill.ppm",10,2);


}
