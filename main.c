#define SDL_MAIN_USE_CALLBACKS 1
#include <SDL3/SDL.h>
#include <SDL3/SDL_main.h>
 
static SDL_Window *window     = NULL;
static SDL_Renderer *renderer = NULL;
 
#define WINDOW_WIDTH  1280
#define WINDOW_HEIGHT 720
#define HALF_WINDOW_WIDTH  (WINDOW_WIDTH  / 2)
#define HALF_WINDOW_HEIGHT (WINDOW_HEIGHT / 2)

typedef struct { float x, y; } Vector2;

typedef struct {
  Vector2 position;
  Vector2 velocity;
  SDL_Color color;
  float lifetime;
  float radius;
  float mass;
} Particle;

#define MAX_PARTICLES 100
Particle particles[MAX_PARTICLES];

SDL_AppResult SDL_AppInit(void **appstate, int argc, char *argv[]) {
  // SDL_SetHint(SDL_HINT_MAIN_CALLBACK_RATE, "120");
  SDL_SetHint(SDL_HINT_RENDER_VSYNC, "1");

  SDL_SetAppMetadata("2D Physics Engine", "0.1", "com.douglasselias.2d-physics-engine");

  if (!SDL_Init(SDL_INIT_VIDEO)) {
    SDL_Log("Couldn't initialize SDL: %s", SDL_GetError());
    return SDL_APP_FAILURE;
  }

  if (!SDL_CreateWindowAndRenderer("examples/renderer/rectangles", WINDOW_WIDTH, WINDOW_HEIGHT, 0, &window, &renderer)) {
    SDL_Log("Couldn't create window/renderer: %s", SDL_GetError());
    return SDL_APP_FAILURE;
  }

  return SDL_APP_CONTINUE;
}

Vector2 spring_end = {HALF_WINDOW_WIDTH, 300};
bool is_dragging = false;


typedef struct {
  float radius;
  Vector2 position;
  Vector2 velocity;
  Vector2 angular_velocity;
} Circle;

Circle a = {150, {HALF_WINDOW_WIDTH, HALF_WINDOW_HEIGHT}, {0,0}};
Circle b = {150, {HALF_WINDOW_WIDTH - 100, -150}, {0,0}};

 /* This function runs when a new event (mouse input, keypresses, etc) occurs. */
SDL_AppResult SDL_AppEvent(void *appstate, SDL_Event *event) {
  if(event->type == SDL_EVENT_QUIT) {
    return SDL_APP_SUCCESS;
  }

  if(event->type == SDL_EVENT_MOUSE_BUTTON_DOWN) {
    is_dragging = true;
    float x, y;
    SDL_GetMouseState(&x, &y);
    spring_end.x = x;
    spring_end.y = y;
  
    // SDL_srand(0);
    // for(int i = 0; i < MAX_PARTICLES; i++) {
    //   particles[i].radius = SDL_randf() * 15;
    //   particles[i].mass = particles[i].radius / 2;
    //   particles[i].lifetime = 5;
    //   particles[i].velocity.x = SDL_randf()*50 * (SDL_randf() > 0.5 ? 1 : -1);
    //   particles[i].velocity.y = SDL_randf()*50 * (SDL_randf() > 0.5 ? 1 : -1);
    //   particles[i].position.x = x;
    //   particles[i].position.y = y;
    //   particles[i].color = (SDL_Color){SDL_randf()*255, SDL_randf()*255, SDL_randf()*255, SDL_ALPHA_OPAQUE};
    // }
  } else {
    is_dragging = false;
  }

  if(event->type == SDL_EVENT_KEY_DOWN) {
    SDL_Scancode key_code = event->key.scancode;
    if(key_code == SDL_SCANCODE_ESCAPE || key_code == SDL_SCANCODE_Q) {
      return SDL_APP_SUCCESS;
    }
  }

  return SDL_APP_CONTINUE;
}


float vector_distance(Vector2 v1, Vector2 v2) {
  float result = sqrtf((v1.x - v2.x)*(v1.x - v2.x) + (v1.y - v2.y)*(v1.y - v2.y));
  return result;
}

float vector_length(Vector2 v) {
  float result = sqrtf((v.x*v.x) + (v.y*v.y));
  return result;
}

Vector2 vector_scale(Vector2 v, float scale) {
  Vector2 result = { v.x*scale, v.y*scale };
  return result;
}

Vector2 normalize_vector(Vector2 v) {
  Vector2 result = { 0 };
  float length = sqrtf((v.x*v.x) + (v.y*v.y));

  if (length > 0) {
    float ilength = 1.0f/length;
    result.x = v.x*ilength;
    result.y = v.y*ilength;
  }

  return result;
}

Vector2 multiply_vectors(Vector2 v1, Vector2 v2) {
  Vector2 result = { v1.x * v2.x, v1.y * v2.y };
  return result;
}

Vector2 sum_vectors(Vector2 v1, Vector2 v2) {
  Vector2 result = { v1.x + v2.x, v1.y + v2.y };
  return result;
}

Vector2 subtract_vectors(Vector2 v1, Vector2 v2) {
  Vector2 result = { v1.x - v2.x, v1.y - v2.y };
  return result;
}

Vector2 subtract_value_from_vector(Vector2 v, float sub) {
  Vector2 result = { v.x - sub, v.y - sub };
  return result;
}

void render_line(Vector2 start, Vector2 end) {
  SDL_RenderLine(renderer, start.x, start.y, end.x, end.y);
}

Vector2 screen_to_cartesian(Vector2 screen) {
  Vector2 cartesian = {0};
  cartesian.x =  screen.x - (WINDOW_WIDTH  / 2);
  cartesian.y = -screen.y + (WINDOW_HEIGHT / 2);
  return cartesian;
}

Vector2 cartesian_to_screen(Vector2 cartesian) {
  Vector2 screen = {0};
  screen.x =  cartesian.x + (WINDOW_WIDTH  / 2);
  screen.y = -cartesian.y + (WINDOW_HEIGHT / 2);
  return screen;
}

float cross_product(Vector2 aaa, Vector2 bbb) {
  return (aaa.x * bbb.y) - (aaa.y * bbb.x);
}

void draw_circle(Vector2 center, float radius, SDL_Color c) {
  SDL_SetRenderDrawColor(renderer, c.r, c.g, c.b, SDL_ALPHA_OPAQUE);

  for(int x = -radius; x <= radius; x++) {
    for(int y = -radius; y <= radius; y++) {
      if(x*x + y*y < radius*radius) {
        SDL_RenderPoint(renderer, center.x + x, center.y + y);
      }
    }
  }
}

Vector2 rotate_vector(Vector2 v, float angle) {
  float s = sinf(angle);
  float c = cosf(angle);
  return (Vector2){ v.x * c - v.y * s, v.x * s + v.y * c };
}

Vector2 apply_drag(Vector2 velocity, float drag_coefficient, float delta_time) {
  float speed = vector_length(velocity);
  if (speed <= 0.0001f) return velocity;  // Avoid division by zero

  float drag_magnitude = drag_coefficient * speed * speed;
  Vector2 drag_force = vector_scale(normalize_vector(velocity), -drag_magnitude);

  velocity.x += drag_force.x * delta_time;
  velocity.y += drag_force.y * delta_time;
  return velocity;
}

Vector2 apply_friction(Vector2 velocity, float friction_coefficient, float normal_force, float delta_time) {
  float speed = vector_length(velocity);
  if (speed < 0.0001f) return (Vector2){0, 0};  // Stop when very slow
  
  float friction_magnitude = friction_coefficient * normal_force;
  Vector2 friction_force = vector_scale(normalize_vector(velocity), -friction_magnitude);

  velocity.x += friction_force.x * delta_time;
  velocity.y += friction_force.y * delta_time;

  // Ensure it doesn't reverse due to friction
  if (vector_length(velocity) > speed) return (Vector2){0, 0};

  return velocity;
}

Vector2 apply_spring_force(Vector2 position, Vector2 anchor, float k, float rest_length, float damping, Vector2 velocity) {
  Vector2 displacement = subtract_vectors(position, anchor);
  
  float length = vector_length(displacement);

  // Normalize displacement to get direction
  Vector2 direction = normalize_vector(displacement);

  // Hooke’s Law: F = -k * (length - rest_length)
  float stretch = length - rest_length;
  Vector2 spring_force = vector_scale(direction, -k * stretch);

  // Damping force: F = -b * velocity
  Vector2 damping_force = vector_scale(velocity, -damping);

  // Total force = Spring force + Damping force
  Vector2 total_force = sum_vectors(spring_force, damping_force);

  return total_force;
}


bool check_collision(Circle *circle_a, Circle *circle_b, float *contact_x, float *contact_y) {
  Vector2 diff = subtract_vectors(circle_b->position, circle_a->position);
  float distance_sq = diff.x * diff.x + diff.y * diff.y;
  float radius_sum = circle_a->radius + circle_b->radius;
  
  if (distance_sq < radius_sum * radius_sum) {  // Collision detected
      float distance = sqrtf(distance_sq);
      float nx = diff.x / distance; // Normalized collision normal
      float ny = diff.y / distance;

      *contact_x = circle_a->position.x + nx * circle_a->radius;
      *contact_y = circle_a->position.y + ny * circle_a->radius;

      return true;
  }

  return false;
}

void resolve_collision(Circle *circle_a, Circle *circle_b) {
  Vector2 diff = subtract_vectors(circle_b->position, circle_a->position);
  float distance = sqrtf(diff.x * diff.x + diff.y * diff.y);
  float nx = diff.x / distance;
  float ny = diff.y / distance;

  Vector2 rel_vel = subtract_vectors(circle_b->velocity, circle_a->velocity);
  float dot = (rel_vel.x * nx + rel_vel.y * ny);

  if (dot > 0) return; // Objects are separating

  float restitution = 0.8f;  // Bounciness
  float impulse = (-(1 + restitution) * dot) / 2;

  circle_a->velocity.x -= impulse * nx;
  circle_a->velocity.y -= impulse * ny;
  circle_b->velocity.x += impulse * nx;
  circle_b->velocity.y += impulse * ny;
}


Vector2 spring_end_vel = {0,0};
static Uint64 last_time = 0;
SDL_AppResult SDL_AppIterate(void *appstate) {
  const Uint64 now = SDL_GetTicks();
  const float delta_time = ((float) (now - last_time)) / 1000.0f;  /* seconds since last iteration */

  Vector2 repulsor = {HALF_WINDOW_WIDTH, HALF_WINDOW_HEIGHT + (HALF_WINDOW_HEIGHT / 2)};
  float repulsor_radius = 90;

  Vector2 attractor = {HALF_WINDOW_WIDTH - (HALF_WINDOW_WIDTH * 0.9), HALF_WINDOW_HEIGHT - (HALF_WINDOW_HEIGHT / 3)};
  float attractor_radius = 90;

  Vector2 forces = {0, 0};
  Vector2 gravity = {0, 70};
  Vector2 random_force = {30, -30};

  for(int i = 0; i < MAX_PARTICLES; i++) {
    if(particles[i].lifetime > 0) {
      forces.x += gravity.x;
      forces.y += gravity.y;

      // Vector2 random_force_acc = {random_force.x / particles[i].mass, random_force.y / particles[i].mass};
      // forces.x += random_force_acc.x;
      // forces.y += random_force_acc.y;

      if(vector_distance(particles[i].position, repulsor) < repulsor_radius * 4) {
        Vector2 cartesian_position = screen_to_cartesian(particles[i].position);
        Vector2 cartesian_repulsor = screen_to_cartesian(repulsor);
        Vector2 repulsion_force = subtract_vectors(cartesian_position, cartesian_repulsor);
        repulsion_force = normalize_vector(repulsion_force);
        float strength = 2300;
        repulsion_force = vector_scale(repulsion_force, strength);
        Vector2 repulsion_force_acc = {repulsion_force.x / particles[i].mass, repulsion_force.y / particles[i].mass};
        forces.x += repulsion_force_acc.x;
        forces.y += repulsion_force_acc.y;
      }

      if(vector_distance(particles[i].position, attractor) < attractor_radius * 4) {
        Vector2 cartesian_position = screen_to_cartesian(particles[i].position);
        Vector2 cartesian_attractor = screen_to_cartesian(attractor);
        Vector2 repulsion_force = subtract_vectors(cartesian_attractor, cartesian_position);
        repulsion_force = normalize_vector(repulsion_force);
        float strength = 2300;
        repulsion_force = vector_scale(repulsion_force, strength);
        Vector2 repulsion_force_acc = {repulsion_force.x / particles[i].mass, repulsion_force.y / particles[i].mass};
        forces.x += repulsion_force_acc.x;
        forces.y += repulsion_force_acc.y;
      }


      /// Integration
      particles[i].velocity.x += forces.x * delta_time;
      particles[i].velocity.y += forces.y * delta_time;

      particles[i].position.x += particles[i].velocity.x * delta_time;
      particles[i].position.y += particles[i].velocity.y * delta_time;

      forces.x = 0;
      forces.y = 0;
    }
  }

  Vector2 anchor = {HALF_WINDOW_WIDTH, 0};
  float spring_mass = 10;

  float k = 50;
  float rest_length = 100;
  float damping = 0.5;
  if(is_dragging) {
    Vector2 spring_force = apply_spring_force(spring_end, anchor, k, rest_length, damping, spring_end_vel);
    spring_end_vel.x += (spring_force.x / spring_mass) * delta_time;
    spring_end_vel.y += (spring_force.y / spring_mass) * delta_time;
    spring_end.x += spring_end_vel.x * delta_time;
    spring_end.x += spring_end_vel.y * delta_time;
  }

  b.velocity.y += gravity.y * delta_time;
  b.position.y += b.velocity.y * delta_time;

  last_time = now;

  SDL_SetRenderDrawColor(renderer, 30, 30, 30, SDL_ALPHA_OPAQUE);
  SDL_RenderClear(renderer);

  float x, y;
  SDL_GetMouseState(&x, &y);
  // b.position.x = x;
  // b.position.y = y;
  Vector2 screen_mouse  = {x, y};
  Vector2 screen_center = {HALF_WINDOW_WIDTH, HALF_WINDOW_HEIGHT};

  Vector2 cartesian_screen = screen_to_cartesian(screen_center);
  Vector2 cartesian_mouse  = screen_to_cartesian(screen_mouse);

  // draw_circle(screen_mouse, 30);
  
  for(int i = 0; i < MAX_PARTICLES; i++) {
    if(particles[i].lifetime > 0)
      draw_circle(particles[i].position, particles[i].radius, particles[i].color);
  }

  draw_circle(repulsor, repulsor_radius, (SDL_Color){200, 200, 200});
  draw_circle(attractor, attractor_radius, (SDL_Color){200, 200, 200});

  SDL_Color reddish_color = {200, 100, 100};
  SDL_SetRenderDrawColor(renderer, reddish_color.r, reddish_color.g, reddish_color.b, SDL_ALPHA_OPAQUE);
  SDL_RenderLine(renderer, anchor.x, anchor.y, spring_end.x, spring_end.y);
  draw_circle(anchor, 10, reddish_color);
  draw_circle(spring_end, 20, reddish_color);

  SDL_Color blueish_color = {70, 70, 200};
  SDL_Color redder_color = {230, 70, 70};
  float contact_a;
  float contact_b;
  if(check_collision(&a, &b, &contact_a, &contact_b)) { 
    draw_circle(a.position, a.radius, redder_color);
    draw_circle(b.position, b.radius, redder_color);
    resolve_collision(&a, &b);
    Vector2 impulse = {10, 10};
    b.velocity.x += impulse.x;
    b.velocity.y += impulse.y;
    float cross = cross_product((Vector2){contact_a, contact_b}, impulse);
    b.angular_velocity.x += cross * delta_time;
    b.angular_velocity.y += cross * delta_time;
  } else {
    draw_circle(a.position, a.radius, blueish_color);
    draw_circle(b.position, b.radius, blueish_color);
  }

  SDL_Color yelloish_color = {170, 170, 30};
  SDL_Color greenish_color = {50, 210, 30};
  Vector2 con = {contact_a, contact_b};
  draw_circle(con, 30, yelloish_color);
  // draw_circle(b.position, b.radius, greenish_color);

  SDL_RenderPresent(renderer);

  return SDL_APP_CONTINUE;
}

void SDL_AppQuit(void *appstate, SDL_AppResult result) {}
