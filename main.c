#include <SDL3/SDL.h>
#include <SDL3/SDL_main.h>
#include <float.h>

#include "src/v2.c"

#define WINDOW_WIDTH 1280
#define WINDOW_HEIGHT 720

typedef SDL_FRect RECT;

typedef enum {
  CIRCLE,
  BOX,
} ShapeType;

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

  SDL_Color color;
  SDL_Color default_color;

  ShapeType shape_type;
} Body;

void get_transformed_vertices(Body* a) {
  for(int i = 0; i < 4; i++) {
    a->transformed_vertices[i] = v2_transform(a->vertices[i], a->position, a->rotation);
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

  body.shape_type = BOX;

  return body;
}

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

  body.shape_type = CIRCLE;

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

float generate_random_radius() {
  return SDL_rand(30) + 15;
}

float generate_random_size() {
  return SDL_rand(100) + 30;
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

bool intersect_polygon(Body a, Body b, V2* normal, float* depth) {
  *normal = (V2){0,0};
  *depth = (float)SDL_MAX_SINT64;

  for(int i = 0; i < 4; i++) {
    V2 va = a.transformed_vertices[i];
    V2 vb = a.transformed_vertices[(i + 1) % 4];

    V2 edge = {vb.x - va.x, vb.y - va.y};
    V2 axis = {-edge.y, edge.x}; // a.k.a normal vector
    axis = v2_normalize(axis);

    float min_a, max_a;
    project_vertices(a.transformed_vertices, axis, &min_a, &max_a);
    float min_b, max_b;
    project_vertices(b.transformed_vertices, axis, &min_b, &max_b);

    if(min_a >= max_b || min_b >= max_a) return false;

    float axis_depth = SDL_min(max_b - min_a, max_a - min_b);
    if(axis_depth < *depth) {
      *depth = axis_depth;
      *normal = axis;
    }
  }

  for(int i = 0; i < 4; i++) {
    V2 va = b.transformed_vertices[i];
    V2 vb = b.transformed_vertices[(i + 1) % 4];

    V2 edge = {vb.x - va.x, vb.y - va.y};
    V2 axis = {-edge.y, edge.x}; // a.k.a normal vector
    axis = v2_normalize(axis);

    float min_a, max_a;
    project_vertices(a.transformed_vertices, axis, &min_a, &max_a);
    float min_b, max_b;
    project_vertices(b.transformed_vertices, axis, &min_b, &max_b);

    if(min_a >= max_b || min_b >= max_a) return false;

    float axis_depth = SDL_min(max_b - min_a, max_a - min_b);
    if(axis_depth < *depth) {
      *depth = axis_depth;
      *normal = axis;
    }
  }

  V2 center_a = get_center_polygon(a.transformed_vertices);
  V2 center_b = get_center_polygon(b.transformed_vertices);

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

  *min = v2_dot(p1, direction);
  *max = v2_dot(p2, direction);

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

bool intersect_circle_polygon(Body a, Body b, V2* normal, float* depth) {
  V2 circle_center = a.position;
  float radius = a.radius;

  *normal = (V2){0,0};
  *depth = (float)SDL_MAX_SINT64;

  for(int i = 0; i < 4; i++) {
    V2 va = b.transformed_vertices[i];
    V2 vb = b.transformed_vertices[(i + 1) % 4];

    V2 edge = {vb.x - va.x, vb.y - va.y};
    V2 axis = {-edge.y, edge.x}; // a.k.a normal vector
    axis = v2_normalize(axis);

    float min_a, max_a;
    project_vertices(b.transformed_vertices, axis, &min_a, &max_a);
    float min_b, max_b;
    project_circle(circle_center, radius, axis, &min_b, &max_b);

    if(min_a >= max_b || min_b >= max_a) return false;

    float axis_depth = SDL_min(max_b - min_a, max_a - min_b);
    if(axis_depth < *depth) {
      *depth = axis_depth;
      *normal = axis;
    }
  }

  int closest_index = find_closest_point_on_polygon(circle_center, b.transformed_vertices);
  V2 closest_point = b.transformed_vertices[closest_index];

  V2 axis = {closest_point.x - circle_center.x, closest_point.y - circle_center.y};
  axis = v2_normalize(axis);

  float min_a, max_a;
  project_vertices(b.transformed_vertices, axis, &min_a, &max_a);
  float min_b, max_b;
  project_circle(circle_center, radius, axis, &min_b, &max_b);

  if(min_a >= max_b || min_b >= max_a) return false;

  float axis_depth = SDL_min(max_b - min_a, max_a - min_b);
  if(axis_depth < *depth) {
    *depth = axis_depth;
    *normal = axis;
  }

  V2 center_polygon = get_center_polygon(b.transformed_vertices);

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

  aabb.min.x = (float)FLT_MAX;
  aabb.min.y = (float)FLT_MAX;

  aabb.max.x = (float)-FLT_MAX;
  aabb.max.y = (float)-FLT_MAX;

  for(int i = 0; i < 4; i++) {
    V2 v = b.vertices[i];
    v.x += b.position.x;
    v.y += b.position.y;

    if(v.x < aabb.min.x) aabb.min.x = v.x;
    if(v.y < aabb.min.y) aabb.min.y = v.y;

    if(v.x > aabb.max.x) aabb.max.x = v.x;
    if(v.y > aabb.max.y) aabb.max.y = v.y;
  }

  return aabb;
}

AABB get_aabb(Body b) {
  switch(b.shape_type) {
    case CIRCLE: return get_aabb_circle(b);
    case BOX:    return get_aabb_polygon(b);
    default:     return (AABB){0};
  }
}

bool intersect_aabb(AABB a, AABB b) {
  if(a.max.x <= b.min.x
  || b.max.x <= a.min.x
  || a.max.y <= b.min.y
  || b.max.y <= a.min.y) return false;

  return true;
}

typedef struct {
  Body* body_a;
  Body* body_b;
  V2 normal;
  float depth;
  V2* contact1;
  V2* contact2;
  int* contact_count;
} Manifold;

void find_contact_point_circles(Body a, Body b, V2* contact_point) {
  V2 direction = v2_sub(b.position, a.position);
  V2 normalized_direction = v2_normalize(direction);
  *contact_point = (V2){
    a.position.x + normalized_direction.x * a.radius,
    a.position.y + normalized_direction.y * a.radius,
  };
}

void point_segment_distance(V2 p, V2 a, V2 b, float* distance_squared, V2* contact_point) {
  V2 ab = v2_sub(b, a);
  V2 ap = v2_sub(p, a);

  float projection = v2_dot(ap, ab);
  float ab_length_squared = v2_length_squared(ab);
  float d = projection / ab_length_squared;

  if(d <= 0) *contact_point = a;
  else if(d >= 1) *contact_point = b;
  else *contact_point = (V2){a.x + ab.x * d, a.y + ab.y * d};

  *distance_squared = v2_distance_squared(p, *contact_point);
}

void find_contact_point_circle_box(Body circle, Body box, V2* contact_point) {
  *contact_point = (V2){0,0};

  float min_distance_squared = FLT_MAX;

  for(int i = 0; i < 4; i++) {
    V2 va = box.transformed_vertices[i];
    V2 vb = box.transformed_vertices[(i + 1) % 4];

    float distance_squared;
    V2 contact;
    point_segment_distance(circle.position, va, vb, &distance_squared, &contact);
    if(distance_squared < min_distance_squared) {
      min_distance_squared = distance_squared;
      *contact_point = contact;
    }
  }
}

bool float_equals(float a, float b) {
  return SDL_fabsf(a - b) < 0.0005f;
  // return SDL_fabsf(a - b) < SDL_FLT_EPSILON;
}

bool v2_equals(V2 a, V2 b) {
  bool x = float_equals(a.x, b.x);
  bool y = float_equals(a.y, b.y);
  return x && y;
}

void find_contact_point_box_box(Body a, Body b, V2* contact1, V2* contact2, int* contact_count) {
  *contact1 = (V2){0,0};
  *contact2 = (V2){0,0};
  *contact_count = 0;

  float min_distance_squared = FLT_MAX;

  for(int i = 0; i < 4; i++) {
    V2 p = a.transformed_vertices[i];

    for(int j = 0; j < 4; j++) {
      V2 va = b.transformed_vertices[j];
      V2 vb = b.transformed_vertices[(j + 1) % 4];

      float distance_squared;
      V2 contact_point;
      point_segment_distance(p, va, vb, &distance_squared, &contact_point);

      if(float_equals(distance_squared, min_distance_squared)) {
        if(!v2_equals(contact_point, *contact1)
        && !v2_equals(contact_point, *contact2)) {
          *contact2 = contact_point;
          *contact_count = 2;
        }
      } else if(distance_squared < min_distance_squared) {
        min_distance_squared = distance_squared;
        *contact_count = 1;
        *contact1 = contact_point;
      }
    }
  }


  for(int i = 0; i < 4; i++) {
    V2 p = b.transformed_vertices[i];

    for(int j = 0; j < 4; j++) {
      V2 va = a.transformed_vertices[j];
      V2 vb = a.transformed_vertices[(j + 1) % 4];

      float distance_squared;
      V2 contact_point;
      point_segment_distance(p, va, vb, &distance_squared, &contact_point);

      if(float_equals(distance_squared, min_distance_squared)) {
        if(!v2_equals(contact_point, *contact1)
        && !v2_equals(contact_point, *contact2)) {
          *contact2 = contact_point;
          *contact_count = 2;
        }
      } else if(distance_squared < min_distance_squared) {
        min_distance_squared = distance_squared;
        *contact_count = 1;
        *contact1 = contact_point;
      }
    }
  }
}

V2 contact_points[1000] = {0};
int contact_index_temp = 0;

void find_contact_points(Body a, Body b, V2* contact1, V2* contact2, int* contact_count) {
  *contact1 = (V2){0,0};
  *contact2 = (V2){0,0};
  *contact_count = 0;

  ShapeType type_a = a.shape_type;
  ShapeType type_b = b.shape_type;

  if(type_a == type_b) {
    switch(type_a) {
      case CIRCLE:
        find_contact_point_circles(a, b, contact1);
        *contact_count = 1;
        break;
      case BOX:
        find_contact_point_box_box(a, b, contact1, contact2, contact_count);
        break;
    }
  } else if(type_a != type_b) {
    switch(type_a) {
      case CIRCLE:
        find_contact_point_circle_box(a, b, contact1);
        *contact_count = 1;
        break;
      case BOX:
        find_contact_point_circle_box(b, a, contact1);
        *contact_count = 1;
        break;
    }
  }

  contact_points[contact_index_temp] = *contact1;
  contact_index_temp = (contact_index_temp+1) % 1000;

  if(!v2_equals(*contact2, (V2){0,0})) {
    contact_points[contact_index_temp] = *contact2;
    contact_index_temp = (contact_index_temp+1) % 1000;
  }
}

int main(int argc, char *argv[]) {
  SDL_Init(SDL_INIT_VIDEO);
  SDL_Window *window = SDL_CreateWindow("2d Physics Engine!", WINDOW_WIDTH, WINDOW_HEIGHT, 0);
  SDL_Renderer *renderer = SDL_CreateRenderer(window, NULL);

  SDL_srand(0);

  bool running = true;
  Uint64 last_time = SDL_GetTicks();

  #define BODIES_COUNT 250
  Body bodies[BODIES_COUNT] = {0};
  int bodies_insert_index = 0;

  #define CONTACTS_COUNT 1000
  Manifold contacts[CONTACTS_COUNT] = {0};
  int contacts_insert_index = 0;

  {
    float density = 1;
    float restitution = 0.5f;
    float width  = WINDOW_WIDTH - 150;
    float height = 100;
    bool is_static = true;
    SDL_Color color = {33, 33, 33, 255};
    bodies[bodies_insert_index] = create_box((V2){WINDOW_WIDTH / 2, WINDOW_HEIGHT - height - 10}, width, height, density, restitution, color, is_static);
    bodies_insert_index++;
  }

  while(running) {
    SDL_Event event;
    while(SDL_PollEvent(&event)) {
      if(event.type == SDL_EVENT_QUIT) running = false;
      if(event.type == SDL_EVENT_MOUSE_BUTTON_UP) {
        bool is_left_click  = event.button.button == SDL_BUTTON_LEFT;
        bool is_right_click = event.button.button == SDL_BUTTON_RIGHT;

        SDL_Color color = random_color_from_raylib_palette();
        float density = 1;
        float restitution = 0.9f;
        float radius = generate_random_radius();
        bool is_static = radius < 20;
        float width  = generate_random_size();
        float height = generate_random_size();

        float x, y;
        SDL_GetMouseState(&x, &y);

        if(is_left_click) {
          bodies[bodies_insert_index] = create_circle((V2){x, y}, radius, density, restitution, color, is_static);
          bodies_insert_index = (bodies_insert_index + 1) % BODIES_COUNT;
        } else if(is_right_click) {
          bodies[bodies_insert_index] = create_box((V2){x, y}, width, height, density, restitution, color, is_static);
          bodies_insert_index = (bodies_insert_index + 1) % BODIES_COUNT;
        }
      }
    }

    Uint64 current_time = SDL_GetTicks();
    float delta_time = (current_time - last_time) / 1000.0f; // Seconds
    last_time = current_time;

    const bool *state = SDL_GetKeyboardState(NULL);
    if(state[SDL_SCANCODE_ESCAPE]) running = false;

    int iterations = 0;
    int total_iterations = 40;
    physics_pass:

    delta_time /= total_iterations;

    //////////////////// Update ///////////////////////

    for(int i = 0; i < BODIES_COUNT; i++) {
      Body* b = &bodies[i];
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
    contacts_insert_index = 0;
    int warning_count = 1;

    for(int i = 0; i < BODIES_COUNT; i++) {
      Body *a = &bodies[i];
      if(a->shape_type == BOX)
        get_transformed_vertices(a);
    }

    for(int i = 0; i < BODIES_COUNT - 1; i++) {
      Body* a = &bodies[i];
      AABB body_a_aabb = get_aabb(*a);

      for(int j = i + 1; j < BODIES_COUNT; j++) {
        Body* b = &bodies[j];
        AABB body_b_aabb = get_aabb(*b);

        if(a->is_static && b->is_static) continue;

        if(!intersect_aabb(body_a_aabb, body_b_aabb)) continue;

        bool collided = false;
        V2 normal = {0,0};
        float depth = 0;

        ShapeType type_a = a->shape_type;
        ShapeType type_b = b->shape_type;

        if(type_a == type_b) {
          switch(type_a) {
            case CIRCLE:
              collided = intersect_circles(*a, *b, &normal, &depth);
              break;
            case BOX:
              collided = intersect_polygon(*a, *b, &normal, &depth);
              break;
          }
        } else if(type_a != type_b) {
          switch(type_a) {
            case CIRCLE:
              collided = intersect_circle_polygon(*a, *b, &normal, &depth);
              break;
            case BOX:
              collided = intersect_circle_polygon(*b, *a, &normal, &depth);
              normal.x = -normal.x; // Flip normal for b->a order
              normal.y = -normal.y;
              break;
          }
        }

        if(collided) {
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

          V2* contact1 = SDL_calloc(sizeof(V2), 1);
          V2* contact2 = SDL_calloc(sizeof(V2), 1);
          int* contact_count = SDL_calloc(sizeof(int), 1);
          find_contact_points(*a, *b, contact1, contact2, contact_count);
          Manifold manifold = {a, b, normal, depth, contact1, contact2, contact_count};
          if(contacts_insert_index + 1 < CONTACTS_COUNT) {
            contacts[contacts_insert_index] = manifold;
            contacts_insert_index++;
          } else {
            SDL_Log("Warning: Not enough size to store all contacts. Warning count: %d\n", warning_count++);
          }
        }
      }
    }

    warning_count = 0;

    for(int i = 0; i < contacts_insert_index; i++) {
      Manifold contact = contacts[i];
      resolve_collision(contact.body_a, contact.body_b, contact.normal, contact.depth);
      SDL_free(contact.contact1);
      SDL_free(contact.contact2);
      SDL_free(contact.contact_count);
    }

    iterations++;
    if(iterations < total_iterations) goto physics_pass;

    /// Soft delete ///
    for(int i = 0; i < BODIES_COUNT; i++) {
      Body* a = &bodies[i];
      if(a->position.x < 0 || a->position.x > WINDOW_WIDTH
      || a->position.y < 0 || a->position.y > WINDOW_HEIGHT) {
        a->position = (V2){-100,-100};
        a->linear_velocity = (V2){0,0};
        a->is_static = true;
      }
    }

    ///////////////// Renderer /////////////////////

    SDL_SetRenderTarget(renderer, NULL);
    SDL_SetRenderDrawColor(renderer, 12, 12, 12, 255);
    SDL_RenderClear(renderer);

    SDL_SetRenderDrawColor(renderer, 255, 0, 0, 255);

    for(int i = 0; i < BODIES_COUNT; i++) {
      if(bodies[i].shape_type == CIRCLE) {
        draw_circle(renderer, bodies[i].position, bodies[i].radius, bodies[i].color);
        bodies[i].color = bodies[i].default_color;
      } else if(bodies[i].shape_type == BOX) {
        SDL_Color c = bodies[i].color;
        SDL_SetRenderDrawColor(renderer, c.r, c.g, c.b, 255);
  
        for(int e = 0; e < 4 - 1; e++) {
          V2 a = bodies[i].transformed_vertices[e];
          for(int w = e + 1; w < 4; w++) {
            V2 b = bodies[i].transformed_vertices[w];
            SDL_RenderLine(renderer, a.x, a.y, b.x, b.y);
          }
        }
  
        bodies[i].color = bodies[i].default_color;
      }
    }

    for(int i = 0; i < 1000; i++) {
      V2 c = contact_points[i];
      draw_circle(renderer, c, 7, (SDL_Color){255,255,255,255});
      contact_points[i] = (V2){-100,-100};
    }
    contact_index_temp = 0;

    SDL_RenderPresent(renderer);
  }

  return 0;
}