extern crate piston_window;

use piston_window::*;

struct Particle {
    position: [f64; 2],
    velocity: [f64; 2],
    force: [f64; 2],
    mass: f64,
    is_draggable: bool,
}

impl Particle {
    fn new(position: [f64; 2], velocity: [f64; 2], force: [f64; 2], mass: f64, is_draggable: bool) -> Self {
        Particle {
            position,
            velocity,
            force,
            mass,
            is_draggable,
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
            damping,
        }
    }
}

struct SoftBody {
    particles: Vec<Particle>,
    springs: Vec<Spring>,
    selected_particle: Option<usize>,
    previous_mouse_pos: Option<[f64; 2]>,
}

impl SoftBody {
    fn new() -> Self {
        SoftBody {
            particles: Vec::new(),
            springs: Vec::new(),
            selected_particle: None,
            previous_mouse_pos: None,
        }
    }

    fn add_particle(&mut self, particle: Particle) {
        self.particles.push(particle);
    }

    fn add_spring(&mut self, spring: Spring) {
        self.springs.push(spring);
    }

    fn is_inside(&self, other_point: &Particle) -> bool {
        let points = &self.particles;

        let mut count = 0;
        for i in 0..points.len() {
            let current_point = &points[i];
            let next_point = &points[(i + 1) % points.len()];

            // Check if the line segment intersects the boundary
            if (current_point.position[1] > other_point.position[1]) != (next_point.position[1] > other_point.position[1]) &&
                other_point.position[0] < (next_point.position[0] - current_point.position[0]) * (other_point.position[1] - current_point.position[1]) / (next_point.position[1] - current_point.position[1]) + current_point.position[0]
            {
                count += 1;
            }
        }

        // If the count is odd, the point is inside the boundary
        count % 2 == 1
    }

    fn update(&mut self, dt: f64, window_size: [f64; 2]) {
        let springs = &self.springs;

        let mut updated_velocities = vec![[0.0, 0.0]; self.particles.len()];

        for (i, particle) in self.particles.iter_mut().enumerate() {
            particle.velocity[0] += (updated_velocities[i][0] + particle.force[0] / particle.mass) * dt;
            particle.velocity[1] += (updated_velocities[i][1] + particle.force[1] / particle.mass) * dt;
        }

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

            let radius = 5.0;

            if particle.position[0] - radius < 0.0 {
                particle.position[0] = radius;
                particle.velocity[0] = -particle.velocity[0];
            }

            if particle.position[0] + radius > window_size[0] {
                particle.position[0] = window_size[0] - radius;
                particle.velocity[0] = -particle.velocity[0];
            }

            if particle.position[1] - radius < 0.0 {
                particle.position[1] = radius;
                particle.velocity[1] = -particle.velocity[1];
            }

            if particle.position[1] + radius > window_size[1] {
                particle.position[1] = window_size[1] - radius;
                particle.velocity[1] = -particle.velocity[1];
            }
        }

        for particle in self.particles.iter_mut() {
            particle.force = [0.0, 0.0];
        }
    }

    fn apply_repulsion(&mut self, other: &Self, repulsion_strength: f64) {
        for particle in &mut self.particles {
            if other.is_inside(particle) {
                // Compute the direction towards the center of the other soft body
                let other_center = other.compute_center();
                let direction = [
                    other_center[0] - particle.position[0],
                    other_center[1] - particle.position[1],
                ];

                // Compute the distance and normalize the direction
                let distance = (direction[0].powi(2) + direction[1].powi(2)).sqrt();
                let direction = [direction[0] / distance, direction[1] / distance];

                // Apply the repulsion force
                particle.velocity[0] -= direction[0] * repulsion_strength;
                particle.velocity[1] -= direction[1] * repulsion_strength;
            }
        }
    }

    fn compute_center(&self) -> [f64; 2] {
        let mut center = [0.0, 0.0];
        for particle in &self.particles {
            center[0] += particle.position[0];
            center[1] += particle.position[1];
        }

        center[0] /= self.particles.len() as f64;
        center[1] /= self.particles.len() as f64;
        center
    }


    fn handle_event(&mut self, e: &Event) {
        if let Some(pos) = e.mouse_cursor_args() {
            self.previous_mouse_pos = Some(pos);
        }

        if let Some(Button::Mouse(MouseButton::Left)) = e.press_args() {
            if let Some(pos) = self.previous_mouse_pos {
                self.selected_particle = self
                    .particles
                    .iter()
                    .position(|particle| {
                        let dx = pos[0] - particle.position[0];
                        let dy = pos[1] - particle.position[1];
                        let distance = (dx * dx + dy * dy).sqrt();
                        distance <= 5.0 && particle.is_draggable
                    });
            }
        }

        if let Some(Button::Mouse(MouseButton::Left)) = e.release_args() {
            self.selected_particle = None;
            self.previous_mouse_pos = None;
        }

        if let Some(pos) = self.previous_mouse_pos {
            if let Some(selected_particle) = self.selected_particle {
                if let Some(particle) = self.particles.get_mut(selected_particle) {
                    particle.position = pos;
                }
            }
        }
    }

    fn render(&self, g: &mut G2d, c: Context) {
        let wireframe_mode = true;
        if wireframe_mode {
            let mut points = Vec::new();

            for particle in &self.particles {
                points.push([particle.position[0], particle.position[1]])
            }

            polygon([1.0, 0.0, 0.0, 1.0], &points, c.transform, g);
        } else {
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
    }
}

fn main() {
    let mut window: PistonWindow =
        WindowSettings::new("Jelly Physics", [800, 600])
            .exit_on_esc(true)
            .build()
            .unwrap();

    let mut soft_body1 = SoftBody::new();
    let mut soft_body2 = SoftBody::new();

    soft_body1.add_particle(Particle::new([200.0, 200.0], [0.0, 0.0], [0.0, 0.0], 1.0, true));
    soft_body1.add_particle(Particle::new([250.0, 200.0], [0.0, 0.0], [0.0, 0.0], 1.0, true));
    soft_body1.add_particle(Particle::new([250.0, 250.0], [0.0, 0.0], [0.0, 0.0], 1.0, true));
    soft_body1.add_particle(Particle::new([200.0, 250.0], [0.0, 0.0], [0.0, 0.0], 1.0, true));

    soft_body1.add_spring(Spring::new(0, 1, 50.0, 3.5, 3.0));
    soft_body1.add_spring(Spring::new(1, 2, 50.0, 3.5, 3.0));
    soft_body1.add_spring(Spring::new(2, 3, 50.0, 3.5, 3.0));
    soft_body1.add_spring(Spring::new(3, 0, 50.0, 3.5, 3.0));
    soft_body1.add_spring(Spring::new(0, 2, 70.7106781187, 8.5, 3.0));
    soft_body1.add_spring(Spring::new(1, 3, 70.7106781187, 8.5, 3.0));

    soft_body2.add_particle(Particle::new([400.0, 200.0], [0.0, 0.0], [0.0, 0.0], 1.0, true));
    soft_body2.add_particle(Particle::new([450.0, 200.0], [0.0, 0.0], [0.0, 0.0], 1.0, true));
    soft_body2.add_particle(Particle::new([450.0, 250.0], [0.0, 0.0], [0.0, 0.0], 1.0, true));
    soft_body2.add_particle(Particle::new([400.0, 250.0], [0.0, 0.0], [0.0, 0.0], 1.0, true));

    soft_body2.add_spring(Spring::new(0, 1, 50.0, 3.5, 3.0));
    soft_body2.add_spring(Spring::new(1, 2, 50.0, 3.5, 3.0));
    soft_body2.add_spring(Spring::new(2, 3, 50.0, 3.5, 3.0));
    soft_body2.add_spring(Spring::new(3, 0, 50.0, 3.5, 3.0));
    soft_body2.add_spring(Spring::new(0, 2, 70.7106781187, 8.5, 3.0));
    soft_body2.add_spring(Spring::new(1, 3, 70.7106781187, 8.5, 3.0));

    while let Some(e) = window.next() {
        if let Some(_) = e.update_args() {
            let deltat_time = 0.01;
            let window_size = [window.size().width as f64, window.size().height as f64];

            let repulsion_constant = 10.0;
            soft_body1.apply_repulsion(&mut soft_body2, repulsion_constant);
            soft_body2.apply_repulsion(&mut soft_body1, repulsion_constant);

            soft_body1.update(deltat_time, window_size);
            soft_body2.update(deltat_time, window_size);
        }

        soft_body1.handle_event(&e);
        soft_body2.handle_event(&e);

        window.draw_2d(&e, |c, g, device| {
            clear([0.5, 0.5, 0.5, 1.0], g);

            soft_body1.render(g, c);
            soft_body2.render(g, c);

            for particle in &soft_body1.particles {
                if soft_body2.is_inside(particle) {
                    ellipse(
                        [1.0, 1.0, 1.0, 1.0], // white color
                        [particle.position[0] - 5.0, particle.position[1] - 5.0, 10.0, 10.0], // circle position and size
                        c.transform,
                        g,
                    );
                }
            }
        });
    }
}
