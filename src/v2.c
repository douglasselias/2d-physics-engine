typedef SDL_FPoint V2;

V2 v2_add(V2 a, V2 b) {
  return (V2){a.x + b.x, a.y + b.y};
}

V2 v2_sub(V2 a, V2 b) {
  return (V2){a.x + b.x, a.y + b.y};
}

float v2_dot(V2 a, V2 b) {
  return a.x * b.x + a.y * b.y;
}

float v2_length(V2 v) {
  return sqrtf(v.x * v.x + v.y * v.y);
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

float v2_distance(V2 a, V2 b) {
  float dx = a.x - b.x;
  float dy = a.y - b.y;
  return sqrt(dx * dx + dy * dy);
}

V2 v2_transform(V2 vertex_target, V2 vertex_base, float angle) {
  float sine   = sin(angle);
  float cosine = cos(angle);

  V2 direction = {
    cosine * vertex_target.x - sine   * vertex_target.y,
    sine   * vertex_target.x + cosine * vertex_target.y,
  };

  return v2_add(vertex_base, direction);
}