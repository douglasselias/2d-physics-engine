#include <SDL3/SDL.h>
#include <SDL3/SDL_main.h>

typedef SDL_FPoint V2;
typedef SDL_FRect RECT;

// typedef enum {
//   CIRCLE,
//   BOX,
// } ShapeType;

typedef struct {
  V2 position;
  V2 linear_velocity;
  float rotation;
  float rotation_velocity;

  float density;
  float mass;
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

  // ShapeType shape_type;

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

void triangulate_box(int triangles[6]) {
  triangles[0] = 0;
  triangles[1] = 1;
  triangles[2] = 2;
  triangles[3] = 0;
  triangles[4] = 2;
  triangles[5] = 3;
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

Body create_box(V2 position, float width, float height, float density, float restitution, SDL_Color color) {
  Body box = {0};
  box.position = position;
  box.density = density;
  box.restitution = restitution;
  box.width = width;
  box.height = height;
  box.area = width * height;
  box.mass = box.area * box.density;
  box.color = color;
  box.default_color = color;
  create_vertices(box.vertices, width, height);
  return box;
}

#define MIN_BODY_SIZE (0.01f * 0.01f)
#define MAX_BODY_SIZE (64.0f * 64.0f)
#define MIN_DENSITY 0.5f // g/cm^3
#define MAX_DENSITY 21.4f

#define WINDOW_WIDTH 1280
#define WINDOW_HEIGHT 720

float calculate_mass(float area, float density) {
  return area * density;
}

float calculate_area_circle(float radius) {
  return 3.14f * radius * radius;
}

Body create_circle(V2 position, float radius, float density, float restitution, SDL_Color color) {
  Body circle = {0};
  circle.position = position;
  circle.radius = radius;
  circle.density = density;
  circle.restitution = restitution;
  circle.area = calculate_area_circle(radius);
  circle.mass = calculate_mass(circle.area, circle.density);
  circle.color = color;
  circle.default_color = color;
  return circle;
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

  *depth /= v2_length(*normal);
  *normal = v2_normalize(*normal);

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
  *min = v2_dot(p1, direction); // added ai
  *max = v2_dot(p2, direction); // added ai

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

    axis = v2_normalize(axis); // added

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

  axis = v2_normalize(axis); // added

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

  // *depth /= v2_length(*normal);
  // *normal = v2_normalize(*normal);

  V2 center_polygon = get_center_polygon(vertices);

  V2 direction = {center_polygon.x - circle_center.x, center_polygon.y - circle_center.y};

  if(v2_dot(direction, *normal) < 0) {
    *normal = (V2){-normal->x, -normal->y};
  }

  return true;
}

int main(int argc, char *argv[]) {
  SDL_Init(SDL_INIT_VIDEO);
  SDL_Window *window = SDL_CreateWindow("2d Physics Engine!", WINDOW_WIDTH, WINDOW_HEIGHT, 0);
  SDL_Renderer *renderer = SDL_CreateRenderer(window, NULL);

  SDL_srand(0);

  bool running = true;
  Uint64 last_time = SDL_GetTicks();

  #define CIRCLES_COUNT 20
  Body circles[CIRCLES_COUNT] = {0};
  for(int i = 0; i < CIRCLES_COUNT; i++) {
    circles[i] = create_circle(generate_random_position(), generate_random_radius(), 0.5, 1, generate_random_color());
  }

  #define RECTS_COUNT 20
  Body rects[RECTS_COUNT] = {0};
  for(int i = 0; i < RECTS_COUNT; i++) {
    rects[i] = create_box(generate_random_position(), generate_random_size(), generate_random_size(), 0.5, 1, generate_random_color());
  }

  while(running) {
    SDL_Event event;
    while(SDL_PollEvent(&event)) {
      if(event.type == SDL_EVENT_QUIT) running = false;
    }

    // Time delta for smooth movement
    Uint64 current_time = SDL_GetTicks();
    float delta_time = (current_time - last_time) / 1000.0f; // Seconds
    last_time = current_time;

    // Keyboard input for square movement
    const bool *state = SDL_GetKeyboardState(NULL);
    if(state[SDL_SCANCODE_ESCAPE]) running = false;

    float speed = 200;
    if(state[SDL_SCANCODE_A]) {
      circles[0].position.x -= speed * delta_time;
      rects[0].position.x   -= speed * delta_time;
    }
    if(state[SDL_SCANCODE_D]) {
      circles[0].position.x += speed * delta_time;
      rects[0].position.x   += speed * delta_time;
    }
    if(state[SDL_SCANCODE_W]) {
      circles[0].position.y -= speed * delta_time;
      rects[0].position.y   -= speed * delta_time;
    }
    if(state[SDL_SCANCODE_S]) {
      circles[0].position.y += speed * delta_time;
      rects[0].position.y   += speed * delta_time;
    }

    for(int i = 0; i < CIRCLES_COUNT - 1; i++) {
      Body* a = &circles[i];
      for(int j = i + 1; j < CIRCLES_COUNT; j++) {
        Body* b = &circles[j];

        V2 normal;
        float depth;
        if(intersect_circles(*a, *b, &normal, &depth)) {
          a->color = (SDL_Color){230, 10, 10, 255};
          b->color = (SDL_Color){230, 10, 10, 255};
          float half_depth = depth / 2;
          a->position.x -= normal.x * half_depth;
          a->position.y -= normal.y * half_depth;
          b->position.x += normal.x * half_depth;
          b->position.y += normal.y * half_depth;
        }
      }
    }

    for(int i = 0; i < RECTS_COUNT; i++) {
      Body *a = &rects[i];
      a->rotation += SDL_PI_F / 2 * delta_time;
      get_transformed_vertices(a);
    }

    for(int i = 0; i < RECTS_COUNT - 1; i++) {
      Body* a = &rects[i];
      for(int j = i + 1; j < RECTS_COUNT; j++) {
        Body* b = &rects[j];

        V2 normal;
        float depth;
        if(intersect_polygon(a->transformed_vertices, b->transformed_vertices, &normal, &depth)) {
          a->color = (SDL_Color){230, 10, 10, 255};
          b->color = (SDL_Color){230, 10, 10, 255};
          float half_depth = depth / 2;
          a->position.x -= normal.x * half_depth;
          a->position.y -= normal.y * half_depth;
          b->position.x += normal.x * half_depth;
          b->position.y += normal.y * half_depth;
        }
      }
    }

    
    for(int i = 0; i < CIRCLES_COUNT; i++) {
      Body* a = &circles[i];
      for(int j = 0; j < RECTS_COUNT; j++) {
        Body* b = &rects[j];

        V2 normal;
        float depth;
        if(intersect_circle_polygon(a->position, a->radius, b->transformed_vertices, &normal, &depth)) {
          a->color = (SDL_Color){230, 10, 10, 255};
          b->color = (SDL_Color){230, 10, 10, 255};
          float half_depth = depth / 2;
          a->position.x -= normal.x * half_depth;
          a->position.y -= normal.y * half_depth;
          b->position.x += normal.x * half_depth;
          b->position.y += normal.y * half_depth;
        }
      }
    }

    ///////////////// Renderer /////////////////////

    SDL_SetRenderTarget(renderer, NULL);
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
    SDL_RenderClear(renderer);

    SDL_SetRenderDrawColor(renderer, 255, 0, 0, 255);
    // RECT rect = {0,0, 100, 200};
    // SDL_RenderFillRect(renderer, &rect);

    for(int i = 0; i < CIRCLES_COUNT; i++) {
      draw_circle(renderer, circles[i].position, circles[i].radius, circles[i].color);
      circles[i].color = circles[i].default_color;
    }

    for(int i = 0; i < RECTS_COUNT; i++) {
      Body r = rects[i];
      // V2 p = r.position;
      // float w = r.width;
      // float h = r.height;
      
      SDL_Color c = r.color;
      SDL_SetRenderDrawColor(renderer, c.r, c.g, c.b, 255);

      // SDL_FRect rect = {p.x, p.y, w, h};
      // SDL_RenderFillRect(renderer, &rect);
      // SDL_RenderLines(renderer,SDL_FPoint *points, count);

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