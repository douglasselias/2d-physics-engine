#include <SDL3/SDL.h>
#include <SDL3/SDL_main.h>

typedef SDL_FPoint V2;
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

  float density;
  float mass;
  float restitution; // value between 0f and 1f
  float area;

  bool is_static;

  float radius;
  float width;
  float height;

  ShapeType shape_type;

  SDL_Color color;
} Body;

// typedef struct {} World;
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

float v2_distance(V2 a, V2 b) {
  float dx = a.x - b.x;
  float dy = a.y - b.y;
  return sqrt(dx * dx + dy * dy);
}

V2 normalize(V2 v) {
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
    *normal = normalize(*normal);
    
    *depth = radii - distance;
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

  #define CIRCLES_COUNT 10
  Body circles[CIRCLES_COUNT] = {0};
  for(int i = 0; i < CIRCLES_COUNT; i++) {
    circles[i] = create_circle(generate_random_position(), generate_random_radius(), 0.5, 1, generate_random_color());
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
    if(state[SDL_SCANCODE_A]) circles[0].position.x -= speed * delta_time;
    if(state[SDL_SCANCODE_D]) circles[0].position.x += speed * delta_time;
    if(state[SDL_SCANCODE_W]) circles[0].position.y -= speed * delta_time;
    if(state[SDL_SCANCODE_S]) circles[0].position.y += speed * delta_time;

    for(int i = 0; i < CIRCLES_COUNT - 1; i++) {
      Body* a = &circles[i];
      for(int j = i + 1; j < CIRCLES_COUNT; j++) {
        Body* b = &circles[j];

        V2 normal;
        float depth;
        if(intersect_circles(*a, *b, &normal, &depth)) {
          float half_depth = depth / 2;
          a->position.x -= normal.x * half_depth;
          a->position.y -= normal.y * half_depth;
          b->position.x += normal.x * half_depth;
          b->position.y += normal.y * half_depth;
        }
      }
    }

    SDL_SetRenderTarget(renderer, NULL);
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
    SDL_RenderClear(renderer);

    SDL_SetRenderDrawColor(renderer, 255, 0, 0, 255);
    // RECT rect = {0,0, 100, 200};
    // SDL_RenderFillRect(renderer, &rect);

    for(int i = 0; i < CIRCLES_COUNT; i++) {
      draw_circle(renderer, circles[i].position, circles[i].radius, circles[i].color);
    }

    SDL_RenderPresent(renderer);
  }

  SDL_DestroyWindow(window);
  SDL_Quit();
  return 0;
}