
#include "field.hpp"
#include "number_theory.hpp"

#include <stdexcept>
#include <utility>

namespace nt = NumberTheory;

template <typename Field> struct EllipticCurve {
  using value_t = typename Field::value_t;
  using field_element_t = FieldElement<Field>;

  struct Point {
    const EllipticCurve &curve;
    field_element_t x, y;
    bool affine;
    Point(const EllipticCurve &curve)
        : curve(curve), x(curve.field(0)), y(curve.field(1)), affine(false) {}
    Point(const EllipticCurve &curve, const field_element_t x,
          const field_element_t y)
        : curve(curve), x(x), y(y), affine(true) {}
    Point(const Point &other) = default;
    Point &operator=(const Point &other) {
      assert(curve == other.curve);
      x = other.x;
      y = other.y;
      affine = other.affine;
      return *this;
    }

    static Point infinity(const EllipticCurve &curve) { return Point(curve); }
    constexpr Point operator-() const {
      return affine ? Point(curve, x, -y) : Point(curve);
    }

    friend Point operator+(const Point &p, const Point &q) {
      if (!p.affine)
        return q;
      if (!q.affine)
        return p;
      if (p.x == q.x && p.y + q.y == 0)
        return infinity(p.curve);
      const field_element_t slope =
          (p.x == q.x) ? (3 * p.x * p.x + p.curve.a) / (2 * p.y)
                       : (q.y - p.y) / (q.x - p.x);
      const field_element_t rx = slope * slope - p.x - q.x;
      const field_element_t ry = slope * (p.x - rx) - p.y;
      return Point(p.curve, rx, ry);
    }

    friend Point operator*(value_t k, Point p) {
      // TODO: Write Montgomery Ladder or w-NAF
      if (k < 0)
        return (-k) * (-p);

      Point result = Point::infinity(p.curve);
      while (k > 0) {
        if (k % 2 == 1)
          result = result + p;
        p = p + p;
        k /= 2;
      }
      return result;
    }
    friend Point operator*(const Point &p, const value_t k) { return k * p; }

    constexpr inline bool operator==(const Point &other) const {
      return x == other.x && y == other.y && affine == other.affine;
    }

    friend std::ostream &operator<<(std::ostream &os, const Point &point) {
      return os << "(" << point.x << " : " << point.y << " : " << point.affine
                << ")";
    }
  };

  using point_t = Point;

  const Field &field;
  const field_element_t a, b;

  EllipticCurve(const Field &field, const field_element_t a,
                const field_element_t b)
      : field(field), a(a), b(b) {
    if (field.characteristic() == 2 || field.characteristic() == 3)
      throw std::runtime_error("EllipticCurve does not currently support "
                               "fields of characteristic 2 or 3");
    const field_element_t discriminant = 4 * a * a * a + 27 * b * b;
    if (discriminant == 0)
      throw std::runtime_error("y^2 = x^3 + " + std::to_string(a.value) +
                               "x + " + std::to_string(b.value) +
                               " defines a singular curve");
  }
  EllipticCurve(const Field &field, const int a, const int b)
      : EllipticCurve(field, field(a), field(b)) {}

  constexpr inline field_element_t f(const field_element_t x) const {
    return x * x * x + a * x + b;
  }
  constexpr inline point_t affine_point(const field_element_t x,
                                        const field_element_t y) const {
    return point_t(*this, x, y);
  }
  constexpr inline point_t affine_point(const value_t x,
                                        const value_t y) const {
    return affine_point(field_element_t(x, field), field_element_t(y, field));
  }
  constexpr inline point_t infinity() const { return point_t(*this); }

  inline std::vector<point_t> points() const {
    std::vector<point_t> result = {Point::infinity(*this)};
    for (value_t x = 0; x < field.p; ++x) {
      const auto x_element = field(x);
      const auto fx = f(x_element);
      if (!nt::is_quadratic_residue(fx.value, field.p))
        continue;

      const auto y_element = field(nt::sqrt_mod(fx.value, field.p));
      const auto p1 = affine_point(x_element, y_element);
      const auto p2 = affine_point(x_element, -y_element);
      result.push_back(p1);
      if (p1 != p2)
        result.push_back(p2);
    }
    return result;
  }

  constexpr inline bool is_point_on_curve(const point_t &point) const {
    if (!point.affine)
      return true;
    return point.y * point.y == f(point.x);
  }

  bool operator==(const EllipticCurve &other) const {
    return field == other.field && a == other.a && b == other.b;
  }

  friend std::ostream &operator<<(std::ostream &os,
                                  const EllipticCurve &curve) {
    return os << "Elliptic curve defined by y^2 = x^3 + " << curve.a.value
              << "*x + " << curve.b.value << " over " << curve.field;
  }
};
