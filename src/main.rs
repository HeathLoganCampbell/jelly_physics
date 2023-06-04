mod square;

extern crate piston_window;

use piston_window::*;
use crate::square::Square;

fn main() {
    let mut window: PistonWindow =
        WindowSettings::new("Hello World!", [512; 2])
            .build().unwrap();

    let square_side_length = 100.0;
    let mut square =
        Square::new(square_side_length, [0.0, 0.0]);

    while let Some(e) = window.next() {
        if let Some(pos) = e.mouse_cursor_args() {
            square.move_to([pos[0], pos[1]]);
        }

        square.render(&mut window, &e);
    }
}