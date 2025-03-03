#include <SDL3/SDL.h>
#include <SDL3/SDL_main.h>

typedef SDL_FPoint V2;
typedef SDL_FRect RECT;

typedef struct {
  V2 position;
  V2 linear_velocity;
  float rotation;
  float rotation_velocity;

  V2 force;

  float density;
  float mass;
  float inv_mass;
  float restitution; // value between 0f and 1f
  float area;

  bool is_static;

  float radius;
  float width;
  float height;

  V2 vertices[4];
  int triangles[6];
  V2 transformed_vertices[4];
  bool transform_update_required;

  SDL_Color color;
  SDL_Color default_color;
} Body;

typedef struct {
  V2 position;
  float sine;
  float cosine;
} Transform;

Transform create_transform(V2 position, float angle) {
  Transform t = {0};
  t.position  = position;
  t.sine      = sin(angle);
  t.cosine    = cos(angle);
  return t;
}

V2 v2_transform(V2 v, Transform t) {
  float rx = t.cosine * v.x - t.sine   * v.y;
  float ry = t.sine   * v.x + t.cosine * v.y;

  float tx = rx + t.position.x;
  float ty = ry + t.position.y;

  return (V2){tx, ty};
}

void get_transformed_vertices(Body* a) {
  Transform t = create_transform(a->position, a->rotation);
  for(int i = 0; i < 4; i++) {
    a->transformed_vertices[i] = v2_transform(a->vertices[i], t);
  }
}

void create_vertices(V2 vertices[4], float width, float height) {
  float left   = -width  / 2;
  float bottom = -height / 2;
  float right  = left   + width;
  float top    = bottom + height;

  vertices[0] = (V2){left,  top};
  vertices[1] = (V2){right, top};
  vertices[2] = (V2){right, bottom};
  vertices[3] = (V2){left,  bottom};
}

Body create_box(V2 position, float width, float height, float density, float restitution, SDL_Color color, bool is_static) {
  Body body = {0};
  body.position = position;
  body.density = density;
  body.restitution = restitution;
  body.width = width;
  body.height = height;
  body.area = width * height;
  // body.mass = body.area * body.density;
  body.mass = 1;
  body.color = color;
  body.default_color = color;
  create_vertices(body.vertices, width, height);

  body.is_static = is_static;
  body.color = color;
  body.default_color = color;

  if(is_static) {
    body.color = (SDL_Color){255, 0, 255, 255};
    body.default_color = (SDL_Color){255, 0, 255, 255};
  }

  body.inv_mass = 1 / body.mass;
  if(is_static) {
    body.inv_mass = 0;
  }

  return body;
}

#define MIN_BODY_SIZE (0.01f * 0.01f)
#define MAX_BODY_SIZE (64.0f * 64.0f)
#define MIN_DENSITY 0.5f // g/cm^3
#define MAX_DENSITY 21.4f

#define WINDOW_WIDTH 1280
#define WINDOW_HEIGHT 720

Body create_circle(V2 position, float radius, float density, float restitution, SDL_Color color, bool is_static) {
  Body body = {0};
  body.position = position;
  body.radius = radius;
  body.density = density;
  body.restitution = restitution;
  body.area = SDL_PI_F * radius * radius;
  // body.mass = body.area * body.density;
  body.mass = 1;

  body.is_static = is_static;
  body.color = color;
  body.default_color = color;

  if(is_static) {
    body.color = (SDL_Color){255, 0, 255, 255};
    body.default_color = (SDL_Color){255, 0, 255, 255};
  }

  body.inv_mass = 1 / body.mass;
  if(is_static) {
    body.inv_mass = 0;
  }

  return body;
}

void draw_circle(SDL_Renderer* renderer, V2 center, float radius, SDL_Color c) {
  SDL_SetRenderDrawColor(renderer, c.r, c.g, c.b, SDL_ALPHA_OPAQUE);

  for(int x = -radius; x <= radius; x++) {
    for(int y = -radius; y <= radius; y++) {
      if(x*x + y*y < radius*radius) {
        SDL_RenderPoint(renderer, center.x + x, center.y + y);
      }
    }
  }
}

SDL_Color generate_random_color() {
  SDL_Color color = {0};
  color.r = (Uint8)SDL_rand(256);
  color.g = (Uint8)SDL_rand(256);
  color.b = (Uint8)SDL_rand(256);
  color.a = 255;
  return color;
}

SDL_Color raylib_palette[] = {
  { 253, 249, 0,   255 }, // YELLOW
  { 255, 203, 0,   255 }, // GOLD
  { 255, 161, 0,   255 }, // ORANGE
  { 230, 41,  55,  255 }, // RED
  { 190, 33,  55,  255 }, // MAROON
  { 0,   228, 48,  255 }, // GREEN
  { 0,   158, 47,  255 }, // LIME
  { 102, 191, 255, 255 }, // SKYBLUE
  { 0,   121, 241, 255 }, // BLUE
  { 0,   82,  172, 255 }, // DARKBLUE
  { 200, 122, 255, 255 }, // PURPLE
  { 135, 60,  190, 255 }, // VIOLET
  { 112, 31,  126, 255 }, // DARKPURPLE
  { 211, 176, 131, 255 }, // BEIGE
  { 127, 106, 79,  255 }, // BROWN
};

SDL_Color random_color_from_raylib_palette() {
  int total_colors = sizeof(raylib_palette) / sizeof(raylib_palette[0]);
  int color_index = SDL_rand(total_colors);
  return raylib_palette[color_index];
}

V2 generate_random_position() {
  V2 position = {0};
  position.x = SDL_rand(WINDOW_WIDTH);
  position.y = SDL_rand(WINDOW_HEIGHT);
  return position;
}

float generate_random_radius() {
  return SDL_rand(30) + 15;
}

float generate_random_size() {
  return SDL_rand(30) + 15;
}

float v2_distance(V2 a, V2 b) {
  float dx = a.x - b.x;
  float dy = a.y - b.y;
  return sqrt(dx * dx + dy * dy);
}

V2 v2_normalize(V2 v) {
  V2 result = {0};
  float length = sqrtf(v.x * v.x + v.y * v.y);

  if(length > SDL_FLT_EPSILON) {
    result.x = v.x / length;
    result.y = v.y / length;
  }

  return result;
}

bool intersect_circles(Body a, Body b, V2* normal, float* depth) {
  float distance = v2_distance(a.position, b.position);
  float radii = a.radius + b.radius;

  *normal = (V2){0,0};
  *depth = 0;

  if(distance >= radii) {
    return false;
  }

  if (distance < SDL_FLT_EPSILON) {
    *normal = (V2){1.0f, 0.0f};
    *depth = radii;
  } else {
    normal->x = b.position.x - a.position.x;
    normal->y = b.position.y - a.position.y;
    *normal = v2_normalize(*normal);
    
    *depth = radii - distance;
  }

  return true;
}

float v2_dot(V2 a, V2 b) {
  return a.x * b.x + a.y * b.y;
}

float v2_length(V2 v) {
  return sqrtf(v.x * v.x + v.y * v.y);
}

V2 get_center_polygon(V2 vertices[4]) {
  V2 center = {0,0};

  for(int i = 0; i < 4; i++) {
    center.x += vertices[i].x;
    center.y += vertices[i].y;
  }

  center.x /= 4;
  center.y /= 4;

  return center;
}

void project_vertices(V2 vertices[4], V2 axis, float* min, float* max) {
   *min = (float)SDL_MAX_SINT64;
   *max = (float)SDL_MIN_SINT64;

   for(int j = 0; j < 4; j++) {
     V2 v = vertices[j];
     float projection = v2_dot(v, axis);

     if(projection < *min) {
       *min = projection;
     }
     
     if(projection > *max) { 
       *max = projection;
     }
   }
}

bool intersect_polygon(V2 vertices_a[4], V2 vertices_b[4], V2* normal, float* depth) {
  *normal = (V2){0,0};
  *depth = (float)SDL_MAX_SINT64;

  for(int i = 0; i < 4; i++) {
    V2 va = vertices_a[i];
    V2 vb = vertices_a[(i + 1) % 4];

    V2 edge = {vb.x - va.x, vb.y - va.y};
    V2 axis = {-edge.y, edge.x}; // a.k.a normal vector
    axis = v2_normalize(axis);

    float min_a, max_a;
    project_vertices(vertices_a, axis, &min_a, &max_a);
    float min_b, max_b;
    project_vertices(vertices_b, axis, &min_b, &max_b);

    if(min_a >= max_b || min_b >= max_a) return false;

    float axis_depth = SDL_min(max_b - min_a, max_a - min_b);
    if(axis_depth < *depth) {
      *depth = axis_depth;
      *normal = axis;
    }
  }

  for(int i = 0; i < 4; i++) {
    V2 va = vertices_b[i];
    V2 vb = vertices_b[(i + 1) % 4];

    V2 edge = {vb.x - va.x, vb.y - va.y};
    V2 axis = {-edge.y, edge.x}; // a.k.a normal vector
    axis = v2_normalize(axis);

    float min_a, max_a;
    project_vertices(vertices_a, axis, &min_a, &max_a);
    float min_b, max_b;
    project_vertices(vertices_b, axis, &min_b, &max_b);

    if(min_a >= max_b || min_b >= max_a) return false;

    float axis_depth = SDL_min(max_b - min_a, max_a - min_b);
    if(axis_depth < *depth) {
      *depth = axis_depth;
      *normal = axis;
    }
  }

  V2 center_a = get_center_polygon(vertices_a);
  V2 center_b = get_center_polygon(vertices_b);

  V2 direction = {center_b.x - center_a.x, center_b.y - center_a.y};

  if(v2_dot(direction, *normal) < 0) {
    *normal = (V2){-normal->x, -normal->y};
  }

  return true;
}

void project_circle(V2 center, float radius, V2 axis, float *min, float* max) {
  V2 direction = v2_normalize(axis);
  V2 direction_and_radius = {direction.x * radius, direction.y * radius};

  V2 p1 = {center.x + direction_and_radius.x, center.y + direction_and_radius.y};
  V2 p2 = {center.x - direction_and_radius.x, center.y - direction_and_radius.y};

  // *min = v2_dot(p1, axis);
  // *max = v2_dot(p2, axis);
  *min = v2_dot(p1, direction); // added, double check if its correct
  *max = v2_dot(p2, direction); // added, double check if its correct

  if(*min > *max) {
    float temp = *min;
    *min = *max;
    *max = temp;
  }
}

int find_closest_point_on_polygon(V2 circle_center, V2 vertices[4]) {
  int index = -1;
  float min_distance = (float)SDL_MAX_SINT64;

  for(int i = 0; i < 4; i++) {
    V2 v = vertices[i];
    float distance = v2_distance(v, circle_center);
    if(distance < min_distance) {
      min_distance = distance;
      index = i;
    }
  }

  return index;
}

bool intersect_circle_polygon(V2 circle_center, float radius, V2 vertices[4], V2*normal, float* depth) {
  *normal = (V2){0,0};
  *depth = (float)SDL_MAX_SINT64;

  for(int i = 0; i < 4; i++) {
    V2 va = vertices[i];
    V2 vb = vertices[(i + 1) % 4];

    V2 edge = {vb.x - va.x, vb.y - va.y};
    V2 axis = {-edge.y, edge.x}; // a.k.a normal vector
    axis = v2_normalize(axis);

    float min_a, max_a;
    project_vertices(vertices, axis, &min_a, &max_a);
    float min_b, max_b;
    project_circle(circle_center, radius, axis, &min_b, &max_b);

    if(min_a >= max_b || min_b >= max_a) return false;

    float axis_depth = SDL_min(max_b - min_a, max_a - min_b);
    if(axis_depth < *depth) {
      *depth = axis_depth;
      *normal = axis;
    }
  }

  int closest_index = find_closest_point_on_polygon(circle_center, vertices);
  V2 closest_point = vertices[closest_index];

  V2 axis = {closest_point.x - circle_center.x, closest_point.y - circle_center.y};
  axis = v2_normalize(axis);

  float min_a, max_a;
  project_vertices(vertices, axis, &min_a, &max_a);
  float min_b, max_b;
  project_circle(circle_center, radius, axis, &min_b, &max_b);

  if(min_a >= max_b || min_b >= max_a) return false;

  float axis_depth = SDL_min(max_b - min_a, max_a - min_b);
  if(axis_depth < *depth) {
    *depth = axis_depth;
    *normal = axis;
  }

  V2 center_polygon = get_center_polygon(vertices);

  V2 direction = {center_polygon.x - circle_center.x, center_polygon.y - circle_center.y};

  if(v2_dot(direction, *normal) < 0) {
    *normal = (V2){-normal->x, -normal->y};
  }

  return true;
}

void resolve_collision(Body* a, Body* b, V2 normal, float depth) {
  V2 relative_velocity = {0,0};

  relative_velocity.x = b->linear_velocity.x - a->linear_velocity.x;
  relative_velocity.y = b->linear_velocity.y - a->linear_velocity.y;

  if(v2_dot(relative_velocity, normal) > 0) return;

  float e = SDL_min(a->restitution, b->restitution);
  float j = -(1 + e) * v2_dot(relative_velocity, normal);

  j /= a->inv_mass + b->inv_mass;

  V2 impulse = {0,0};
  impulse.x = j * normal.x;
  impulse.y = j * normal.y;

  a->linear_velocity.x -= impulse.x * a->inv_mass;
  a->linear_velocity.y -= impulse.y * a->inv_mass;

  b->linear_velocity.x += impulse.x * b->inv_mass;
  b->linear_velocity.y += impulse.y * b->inv_mass;
}

float wrap_value(float value, float min, float max) {
  float result = value - (max - min) * floorf((value - min)/(max - min));
  return result;
}

V2 gravity = {0, 9.81f * 37500};

typedef struct {
  V2 min, max;
} AABB;

AABB get_aabb_circle(Body b) {
  AABB aabb = {0};

  aabb.min.x = b.position.x - b.radius;
  aabb.min.y = b.position.y - b.radius;

  aabb.max.x = b.position.x + b.radius;
  aabb.max.y = b.position.y + b.radius;

  return aabb;
}

AABB get_aabb_polygon(Body b) {
  AABB aabb = {0};

  aabb.min.x = (float)SDL_MAX_SINT64;
  aabb.min.y = (float)SDL_MAX_SINT64;

  aabb.max.x = (float)SDL_MIN_SINT64;
  aabb.max.y = (float)SDL_MIN_SINT64;

  for(int i = 0; i < 4; i++) {
    V2 v = b.vertices[i];

    if(v.x < aabb.min.x) aabb.min.x = v.x;
    if(v.y < aabb.min.y) aabb.min.y = v.y;

    if(v.x > aabb.max.x) aabb.max.x = v.x;
    if(v.y > aabb.max.y) aabb.max.y = v.y;
  }

  return aabb;
}

int main(int argc, char *argv[]) {
  SDL_Init(SDL_INIT_VIDEO);
  SDL_Window *window = SDL_CreateWindow("2d Physics Engine!", WINDOW_WIDTH, WINDOW_HEIGHT, 0);
  SDL_Renderer *renderer = SDL_CreateRenderer(window, NULL);

  SDL_srand(0);

  bool running = true;
  Uint64 last_time = SDL_GetTicks();

  #define CIRCLES_COUNT 50
  Body circles[CIRCLES_COUNT] = {0};
  int circles_insert_index = 0;

  #define RECTS_COUNT 50
  Body rects[RECTS_COUNT] = {0};
  int rects_insert_index = 0;

  {
    float density = 1;
    float restitution = 0.5f;
    float width  = WINDOW_WIDTH - 150;
    float height = 100;
    bool is_static = true;
    SDL_Color color = {33, 33, 33, 255};
    rects[rects_insert_index] = create_box((V2){WINDOW_WIDTH / 2, WINDOW_HEIGHT - height - 10}, width, height, density, restitution, color, is_static);
    rects_insert_index++;
  }

  while(running) {
    SDL_Event event;
    while(SDL_PollEvent(&event)) {
      if(event.type == SDL_EVENT_QUIT) running = false;
      if(event.type == SDL_EVENT_MOUSE_BUTTON_UP) {
        bool is_left_click  = event.button.button == SDL_BUTTON_LEFT;
        bool is_right_click = event.button.button == SDL_BUTTON_RIGHT;

        if(is_left_click) {
          float density = 1;
          float restitution = 0.9f;
          float radius = generate_random_radius();
          bool is_static = radius < 20;
          SDL_Color color = random_color_from_raylib_palette();
          float x, y;
          SDL_GetMouseState(&x, &y);
          circles[circles_insert_index] = create_circle((V2){x, y}, radius, density, restitution, color, is_static);
          circles_insert_index = (circles_insert_index + 1) % CIRCLES_COUNT;
        } else if(is_right_click) {
          float density = 1;
          float restitution = 0.3f;
          float width  = generate_random_size();
          float height = generate_random_size();
          bool is_static = width < 20;
          SDL_Color color = random_color_from_raylib_palette();
          float x, y;
          SDL_GetMouseState(&x, &y);
          rects[rects_insert_index] = create_box((V2){x, y}, width, height, density, restitution, color, is_static);
          rects_insert_index = (rects_insert_index + 1) % RECTS_COUNT;
        }
      }
    }

    Uint64 current_time = SDL_GetTicks();
    float delta_time = (current_time - last_time) / 1000.0f; // Seconds
    // SDL_Log("DT: %.5f\n", delta_time);
    last_time = current_time;

    const bool *state = SDL_GetKeyboardState(NULL);
    if(state[SDL_SCANCODE_ESCAPE]) running = false;

    int iterations = 0;
    int total_iterations = 40;
    physics_pass:

    delta_time /= total_iterations;

    //////////////////// Update ///////////////////////
    for(int i = 0; i < CIRCLES_COUNT; i++) {
      Body* b = &circles[i];
      if(b->is_static) continue;

      // b->force.y += gravity.y;

      // V2 acc = {
      //   b->force.x / b->mass,
      //   b->force.y / b->mass
      // };
      // b->linear_velocity.x += acc.x * delta_time;
      // b->linear_velocity.y += acc.y * delta_time;

      b->linear_velocity.y += gravity.y * delta_time;

      b->position.x += b->linear_velocity.x * delta_time;
      b->position.y += b->linear_velocity.y * delta_time;

      b->rotation += b->rotation_velocity * delta_time;

      /// Don't forget to reset the forces!
      b->force = (V2){0,0};
    }

    for(int i = 0; i < RECTS_COUNT; i++) {
      Body* b = &rects[i];
      if(b->is_static) continue;

      // b->force.y += gravity.y;

      // V2 acc = {
      //   b->force.x / b->mass,
      //   b->force.y / b->mass
      // };

      // b->linear_velocity.x += acc.x * delta_time;
      // b->linear_velocity.y += acc.y * delta_time;

      b->linear_velocity.y += gravity.y * delta_time;

      b->position.x += b->linear_velocity.x * delta_time;
      b->position.y += b->linear_velocity.y * delta_time;

      b->rotation += b->rotation_velocity * delta_time;

      /// Don't forget to reset the forces!
      b->force = (V2){0,0};
    }

    /////////////////////////// Collisions //////////////////////////////
    for(int i = 0; i < CIRCLES_COUNT - 1; i++) {
      Body* a = &circles[i];
      for(int j = i + 1; j < CIRCLES_COUNT; j++) {
        Body* b = &circles[j];

        if(a->is_static && b->is_static) continue;

        V2 normal;
        float depth;
        if(intersect_circles(*a, *b, &normal, &depth)) {
          a->color = (SDL_Color){230, 10, 10, 255};
          b->color = (SDL_Color){230, 10, 10, 255};

          if(a->is_static) {
            b->position.x += normal.x * depth;
            b->position.y += normal.y * depth;
          } else if(b->is_static) {
            a->position.x -= normal.x * depth;
            a->position.y -= normal.y * depth;
          } else {
            float half_depth = depth / 2;
            a->position.x -= normal.x * half_depth;
            a->position.y -= normal.y * half_depth;
            b->position.x += normal.x * half_depth;
            b->position.y += normal.y * half_depth;
          }

          resolve_collision(a, b, normal, depth);
        }
      }
    }

    for(int i = 0; i < RECTS_COUNT; i++) {
      Body *a = &rects[i];
      get_transformed_vertices(a);
    }

    for(int i = 0; i < RECTS_COUNT - 1; i++) {
      Body* a = &rects[i];
      for(int j = i + 1; j < RECTS_COUNT; j++) {
        Body* b = &rects[j];

        if(a->is_static && b->is_static) continue;

        V2 normal;
        float depth;
        if(intersect_polygon(a->transformed_vertices, b->transformed_vertices, &normal, &depth)) {
          a->color = (SDL_Color){230, 10, 10, 255};
          b->color = (SDL_Color){230, 10, 10, 255};

          if(a->is_static) {
            b->position.x += normal.x * depth;
            b->position.y += normal.y * depth;
          } else if(b->is_static) {
            a->position.x -= normal.x * depth;
            a->position.y -= normal.y * depth;
          } else {
            float half_depth = depth / 2;
            a->position.x -= normal.x * half_depth;
            a->position.y -= normal.y * half_depth;
            b->position.x += normal.x * half_depth;
            b->position.y += normal.y * half_depth;
          }
        }
      }
    }

    for(int i = 0; i < CIRCLES_COUNT; i++) {
      Body* a = &circles[i];
      for(int j = 0; j < RECTS_COUNT; j++) {
        Body* b = &rects[j];

        if(a->is_static && b->is_static) continue;

        V2 normal;
        float depth;
        if(intersect_circle_polygon(a->position, a->radius, b->transformed_vertices, &normal, &depth)) {
          a->color = (SDL_Color){230, 10, 10, 255};
          b->color = (SDL_Color){230, 10, 10, 255};

          if(a->is_static) {
            b->position.x += normal.x * depth;
            b->position.y += normal.y * depth;
          } else if(b->is_static) {
            a->position.x -= normal.x * depth;
            a->position.y -= normal.y * depth;
          } else {
            float half_depth = depth / 2;
            a->position.x -= normal.x * half_depth;
            a->position.y -= normal.y * half_depth;
            b->position.x += normal.x * half_depth;
            b->position.y += normal.y * half_depth;
          }

          resolve_collision(a, b, normal, depth);
        }
      }
    }

    iterations++;
    if(iterations < total_iterations) goto physics_pass;

    /// Wrap screen ///
    for(int i = 0; i < CIRCLES_COUNT; i++) {
      Body* a = &circles[i];
      if(a->position.x < 0 || a->position.x > WINDOW_WIDTH
      || a->position.y < 0 || a->position.y > WINDOW_HEIGHT) {
        a->position = (V2){-100,-100};
        a->linear_velocity = (V2){0,0};
        a->is_static = true;
      }
    }

    for(int j = 0; j < RECTS_COUNT; j++) {
      Body* a = &rects[j];
      if(a->position.x < 0 || a->position.x > WINDOW_WIDTH
      || a->position.y < 0 || a->position.y > WINDOW_HEIGHT) {
        a->position = (V2){-100,-100};
        a->linear_velocity = (V2){0,0};
        a->is_static = true;
      }
    }

    ///////////////// Renderer /////////////////////

    SDL_SetRenderTarget(renderer, NULL);
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
    SDL_RenderClear(renderer);

    SDL_SetRenderDrawColor(renderer, 255, 0, 0, 255);

    for(int i = 0; i < CIRCLES_COUNT; i++) {
      draw_circle(renderer, circles[i].position, circles[i].radius, circles[i].color);
      circles[i].color = circles[i].default_color;
    }

    for(int i = 0; i < RECTS_COUNT; i++) {
      Body r = rects[i];
      SDL_Color c = r.color;
      SDL_SetRenderDrawColor(renderer, c.r, c.g, c.b, 255);

      for(int e = 0; e < 4 - 1; e++) {
        V2 a = r.transformed_vertices[e];
        for(int w = e + 1; w < 4; w++) {
          V2 b = r.transformed_vertices[w];
          SDL_RenderLine(renderer, a.x, a.y, b.x, b.y);
        }
      }

      rects[i].color = rects[i].default_color;
    }

    SDL_RenderPresent(renderer);
  }

  SDL_DestroyWindow(window);
  SDL_Quit();
  return 0;
}