extern crate piston_window;

use piston_window::*;

fn main() {
    let mut window: PistonWindow =
        WindowSettings::new("Hello World!", [512; 2])
            .build().unwrap();

    let mut circle_pos: [f64; 2] = [0.0, 0.0];

    while let Some(e) = window.next() {
        if let Some(pos) = e.mouse_cursor_args() {
            circle_pos = pos;
        }

        window.draw_2d(&e, |c, g, _| {
            clear([0.5, 0.5, 0.5, 1.0], g);

            ellipse( [0.0, 0.0, 1.0, 1.0],
                     [circle_pos[0] - 5.0, circle_pos[1] - 5.0, 10.0, 10.0],
                     c.transform,
                     g
            );

            rectangle([1.0, 0.0, 0.0, 1.0], // red
                      [0.0, 0.0, 100.0, 100.0], // rectangle
                      c.transform, g);
        });
    }
}