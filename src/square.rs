use piston_window::*;

pub  struct Point {
    x: f64,
    y: f64,
}

pub struct Square {
    side_length: f64,
    position: [f64; 2],
}

impl Square {
    pub fn new(side_length: f64, position: [f64; 2]) -> Self {
        Square {
            side_length,
            position,
        }
    }

    fn generate_points(&self) -> Vec<Point> {
        let point_a = Point { x: 0.0, y: 0.0 };
        let point_b = Point { x: self.side_length, y: 0.0 };
        let point_c = Point { x: self.side_length, y: self.side_length };
        let point_d = Point { x: 0.0, y: self.side_length };

        vec![point_a, point_b, point_c, point_d]
    }

    pub fn move_to(&mut self, position: [f64; 2]) {
        self.position = position;

        if self.position[1] > 500 as f64 - self.side_length
        {
            self.position[1] = 500 as f64 - self.side_length;
        }
    }

    pub fn render(&self, window: &mut PistonWindow, e: &Event) {
        let square_points = self.generate_points();

        window.draw_2d(e, |c, g, _| {
            clear([0.5, 0.5, 0.5, 1.0], g);

            for point in &square_points {
                ellipse(
                    [1.0, 1.0, 1.0, 1.0], // blue color
                    [point.x + self.position[0], point.y + self.position[1], 10.0, 10.0], // circle position and size
                    c.transform,
                    g,
                );
            }

            for i in 0..square_points.len() {
                let current_point = &square_points[i];
                let next_point = &square_points[(i + 1) % square_points.len()];

                line(
                    [1.0, 1.0, 1.0, 1.0], // red color
                    1.0, // line width
                    [current_point.x + self.position[0] + 5.0, current_point.y + self.position[1] + 5.0, next_point.x + self.position[0] + 5.0, next_point.y + self.position[1] + 5.0], // line start and end points
                    c.transform,
                    g,
                );
            }
        });
    }
}