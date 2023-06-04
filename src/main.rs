extern crate piston_window;

use std::io::Write;
use piston_window::*;

struct Particle {
    position: [f64; 2],
    velocity: [f64; 2],
    mass: f64,
}

impl Particle {
    fn new(position: [f64; 2], velocity: [f64; 2], mass: f64) -> Self {
        Particle {
            position,
            velocity,
            mass,
        }
    }
}

struct Spring {
    particle1_index: usize,
    particle2_index: usize,
    rest_length: f64,
    stiffness: f64,
    damping: f64,
}

impl Spring {
    fn new(particle1_index: usize, particle2_index: usize, rest_length: f64, stiffness: f64, damping: f64) -> Self {
        Spring {
            particle1_index,
            particle2_index,
            rest_length,
            stiffness,
            damping
        }
    }
}

struct SoftBody {
    particles: Vec<Particle>,
    springs: Vec<Spring>,
}

impl SoftBody {
    fn new() -> Self {
        SoftBody {
            particles: Vec::new(),
            springs: Vec::new(),
        }
    }

    fn add_particle(&mut self, particle: Particle) {
        self.particles.push(particle);
    }

    fn add_spring(&mut self, spring: Spring) {
        self.springs.push(spring);
    }

    fn update(&mut self, dt: f64) {
        let springs = &self.springs;

        let mut updated_velocities = vec![[0.0, 0.0]; self.particles.len()];

        for spring in springs {
            let particle1 = &self.particles[spring.particle1_index];
            let particle2 = &self.particles[spring.particle2_index];

            let displacement = [
                particle2.position[0] - particle1.position[0],
                particle2.position[1] - particle1.position[1],
            ];
            let distance = (displacement[0].powi(2) + displacement[1].powi(2)).sqrt();
            let force_magnitude = spring.stiffness * (distance - spring.rest_length);

            let force_direction = [
                displacement[0] / distance,
                displacement[1] / distance,
            ];
            let force = [
                force_direction[0] * force_magnitude,
                force_direction[1] * force_magnitude,
            ];

            updated_velocities[spring.particle1_index][0] += force[0] / particle1.mass;
            updated_velocities[spring.particle1_index][1] += force[1] / particle1.mass;

            updated_velocities[spring.particle2_index][0] -= force[0] / particle2.mass;
            updated_velocities[spring.particle2_index][1] -= force[1] / particle2.mass;
        }

        for (i, particle) in self.particles.iter_mut().enumerate() {
            particle.velocity[0] += updated_velocities[i][0];
            particle.velocity[1] += updated_velocities[i][1];

            // Apply damping force
            particle.velocity[0] -= particle.velocity[0] * springs[i].damping * dt;
            particle.velocity[1] -= particle.velocity[1] * springs[i].damping * dt;

            particle.position[0] += particle.velocity[0] * dt;
            particle.position[1] += particle.velocity[1] * dt;
        }
    }

    fn render(&self, e: &Event, g: &mut G2d, c: Context) {
        for particle in &self.particles {
            ellipse(
                [1.0, 1.0, 1.0, 1.0], // white color
                [particle.position[0] - 5.0, particle.position[1] - 5.0, 10.0, 10.0], // circle position and size
                c.transform,
                g,
            );
        }

        for spring in &self.springs {
            let particle1 = &self.particles[spring.particle1_index];
            let particle2 = &self.particles[spring.particle2_index];

            line(
                [1.0, 1.0, 1.0, 1.0], // white color
                1.0, // line width
                [particle1.position[0], particle1.position[1], particle2.position[0], particle2.position[1]], // line start and end points
                c.transform,
                g,
            );
        }
    }

    fn handle_mouse_drag(&mut self, mouse_x: f64, mouse_y: f64) {
        let closest_particle_index = self
            .particles
            .iter()
            .enumerate()
            .min_by(|(_, p1), (_, p2)| {
                let dx1 = p1.position[0] - mouse_x;
                let dy1 = p1.position[1] - mouse_y;
                let dx2 = p2.position[0] - mouse_x;
                let dy2 = p2.position[1] - mouse_y;
                (dx1 * dx1 + dy1 * dy1)
                    .partial_cmp(&(dx2 * dx2 + dy2 * dy2))
                    .unwrap()
            })
            .map(|(index, _)| index);

        if let Some(index) = closest_particle_index {
            self.particles[index].position = [mouse_x, mouse_y];
        }
    }
}

fn main() {
    let mut window: PistonWindow =
        WindowSettings::new("Soft Body Engine", [800, 600])
            .exit_on_esc(true)
            .build()
            .unwrap();

    let mut soft_body = SoftBody::new();

    // Add particles to the soft body
    soft_body.add_particle(Particle::new([100.0, 100.0], [0.0, 0.0], 5.0));
    soft_body.add_particle(Particle::new([100.0, 200.0], [0.0, 0.0], 5.0));
    soft_body.add_particle(Particle::new([200.0, 100.0], [0.0, 0.0], 5.0));
    soft_body.add_particle(Particle::new([200.0, 200.0], [0.0, 0.0], 5.0));

    // Add springs connecting the particles
    soft_body.add_spring(Spring::new(0, 1, 100.0, 0.7, 3.0));
    soft_body.add_spring(Spring::new(0, 2, 100.0, 0.7, 3.0));
    soft_body.add_spring(Spring::new(1, 3, 100.0, 0.7, 3.0));
    soft_body.add_spring(Spring::new(2, 3, 100.0, 0.7, 3.0));

    soft_body.add_spring(Spring::new(0, 3, 141.421, 0.7, 3.0));
    soft_body.add_spring(Spring::new(2, 1, 141.421, 0.7, 3.0));

    // Square 2
    let mut soft_body2 = SoftBody::new();

    // Add particles to the soft body
    soft_body2.add_particle(Particle::new([100.0, 100.0], [0.0, 0.0], 5.0));
    soft_body2.add_particle(Particle::new([100.0, 200.0], [0.0, 0.0], 5.0));
    soft_body2.add_particle(Particle::new([200.0, 100.0], [0.0, 0.0], 5.0));
    soft_body2.add_particle(Particle::new([200.0, 200.0], [0.0, 0.0], 5.0));

    // Add springs connecting the particles
    soft_body2.add_spring(Spring::new(0, 1, 100.0, 0.7, 3.0));
    soft_body2.add_spring(Spring::new(0, 2, 100.0, 0.7, 3.0));
    soft_body2.add_spring(Spring::new(1, 3, 100.0, 0.7, 3.0));
    soft_body2.add_spring(Spring::new(2, 3, 100.0, 0.7, 3.0));

    soft_body2.add_spring(Spring::new(0, 3, 141.421, 0.7, 3.0));
    soft_body2.add_spring(Spring::new(2, 1, 141.421, 0.7, 3.0));

    let mut dragging_particle = false;

    while let Some(e) = window.next() {
        if let Some(args) = e.mouse_cursor_args() {
            let (mouse_x, mouse_y) = (args[0], args[1]);

            if dragging_particle {
                soft_body.handle_mouse_drag(mouse_x, mouse_y);
            }
        }

        if let Some(Button::Mouse(MouseButton::Left)) = e.press_args() {
            dragging_particle = true;
        }

        if let Some(Button::Mouse(MouseButton::Left)) = e.release_args() {
            dragging_particle = false;
        }

        if let Some(_) = e.update_args() {
            let dt = 0.01;

            soft_body.update(dt);
            soft_body2.update(dt);
        }

        window.draw_2d(&e, |c, g, device| {
            clear([0.5, 0.5, 0.5, 1.0], g);

            soft_body.render(&e, g, c);
            soft_body2.render(&e, g, c);
        });
    }
}
