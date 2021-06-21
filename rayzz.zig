const std = @import("std");

const output_fname = "foo.tga";

// Canvas width x height - canvas is the image we produce, in pixels.
const width: u16 = 640;
const height: u16 = 480;

const Vec3 = struct {
    x: f32 = 0.0,
    y: f32 = 0.0,
    z: f32 = 0.0,

    // Handy shortcut.
    pub fn new(x: f32, y: f32, z: f32) Vec3 {
        return Vec3{ .x = x, .y = y, .z = z };
    }

    // Dot product is just the properties multiplied by each other.
    pub fn dot(self: Vec3, other: Vec3) f32 {
        return self.x * other.x + self.y * other.y + self.z * other.z;
    }

    // The difference of all components of 'other' from 'self'.
    // If the two inputs are points, the output is a vector.
    pub fn minus(self: Vec3, other: Vec3) Vec3 {
        return Vec3.new(
            self.x - other.x,
            self.y - other.y,
            self.z - other.z,
        );
    }

    // Multiply all components by scalar value.
    pub fn scale(self: Vec3, scalar: f32) Vec3 {
        return Vec3.new(self.x * scalar, self.y * scalar, self.z * scalar);
    }
};

const Sphere = struct {
    center: Vec3,
    radius: f32,
    color: Color
};

const Color = struct {
    r: u8 = 100,
    g: u8 = 100,
    b: u8 = 100,
};

// The scene is a list of spheres.
const spheres = [_]Sphere{
    Sphere{
        .center = Vec3.new(0, -0.5, 3),
        .radius = 1,
        .color = Color{ .r = 171, .g = 222, .b = 20 },
    },
    Sphere{
        .center = Vec3.new(0.3, 0.4, 2.5),
        .radius = 0.3,
        .color = Color{ .r = 235, .g = 107, .b = 30 },
    },
    Sphere{
        .center = Vec3.new(-0.6, 0.35, 2.7),
        .radius = 0.4,
        .color = Color{ .r = 245, .g = 188, .b = 10 },
    },
    Sphere{
        .center = Vec3.new(2, 0, 4),
        .radius = 1,
        .color = Color{ .r = 9, .g = 203, .b = 235 },
    },
};

// The "camera" is composed of a viewport and an origin. We trace
// rays from the origin point through points in the viewport into
// the scene.
//
//           /|    ^
//          / |    |
// origin> *  | vp_height
//          \ |    |
//           \|    v
//         <--> vp_dist
//
const vp_width: f32 = 1.0;
const vp_height: f32 = 1.0;
const vp_dist: f32 = 1.0;
const origin_pt = Vec3.new(0, 0, 0);
// Convert from canvas to viewport scale (canvas is based on pixels)
const scale_x: f32 = vp_width / @as(f32, width);
const scale_y: f32 = vp_height / @as(f32, width);

pub fn main() !void {
    //
    // . -x,y |
    //        |          The scene has origin 0,0
    //        | 0,0      and positive and negative
    // -------+------    x and y values above and
    //        |          below.
    //        |
    //        |     . x,-y
    //
    // For every pixel, we plot a point on the viewport, starting
    // at the top left.
    const x_start: i32 = -1 * (@as(i64, width) / 2);
    const y_start: i32 = (@as(i64, height) / 2);
    const x_end = x_start * -1; // negative to positive
    const y_end = y_start * -1; // positive to negative

    var y = y_start;
    while (y > y_end) : (y -= 1) { // down from y to -y

        var x = x_start;
        while (x < x_end) : (x += 1) { // up from -x to x

            // Make a 3D point on the viewport plane corresponding with
            // the canvas pixel we wish to draw.
            const viewport_pt = Vec3{
                .x = @intToFloat(f32, x) * scale_x,
                .y = @intToFloat(f32, y) * scale_y,
                .z = vp_dist,
            };

            // Trace a ray from the origin, through the viewport point and get
            // back a color.
            const color = traceRay(origin_pt, viewport_pt, 1, std.math.f32_max);
            putPixel(x, y, color);
        }
    }

    try writeTga();
}

// Trace a ray!
//
//  o--vp-----------> ???
//
fn traceRay(my_origin_pt: Vec3, viewport_pt: Vec3, min: f32, max: f32) Color {
    var closest_point = max;
    var closest_sphere: ?Sphere = null;

    for (spheres) |sphere| {
        // Get the scalar values needed to "scale" the ray vector so that it
        // intersects with the sphere. The values of t1 and t2 give us the
        // relative distance to the near and far surfaces of the sphere.
        var t1: f32 = undefined;
        var t2: f32 = undefined;
        getRaySphereIntersection(my_origin_pt, viewport_pt, sphere, &t1, &t2);

        // See if we have a new closest point and/or closest sphere.
        if (t1 > min and t1 < max and t1 < closest_point) {
            closest_point = t1;
            closest_sphere = sphere;
        }
        if (t2 > min and t2 < max and t2 < closest_point) {
            closest_point = t2;
            closest_sphere = sphere;
        }
    }

    // Return the color of the closest sphere on this ray, if any.
    if (closest_sphere) |sphere| {
        return sphere.color;
    }

    // No sphere, return background (default) color
    return Color{};
}

// Insersect a sphere!
//
//                "boop"   "boop"
//  o--vp-----------| sphere |---->
//
fn getRaySphereIntersection(
    my_origin_pt: Vec3,
    viewport_pt: Vec3,
    sphere: Sphere,
    t1: *f32,
    t2: *f32,
) void {
    // Note: This algorithm has been algebraically transformed beyond any
    // recognition, but the principle is that we find a solution for these two
    // equations:
    //
    //     |P-C| = r      - The sphere's surface (|P-C| is distance between)
    //     O + t(D) = P   - A ray from the origin, scaled by scalar t
    //
    //      - P is a point on the surface of the sphere
    //      - C is the center of the sphere
    //      - r is the radius of the sphere
    //      - O is the origin of our ray
    //      - D is the vector representing our ray
    //      - t is the scalar value that we're trying to solve for
    //
    // The solution is a quadratic equation and t can be two values: t1 is
    // the intersection at the front surface and t2 is the intersection at the
    // back of the sphere.
    //
    const origin_to_center = my_origin_pt.minus(sphere.center);
    
    // use the expected a,b,c variables for the quadratic equation:
    const a = viewport_pt.dot(viewport_pt);
    const b = origin_to_center.dot(viewport_pt) * 2;
    const c = origin_to_center.dot(origin_to_center) - sphere.radius * sphere.radius;

    // You can determine the "roots" of a quadratic equation by examining the
    // discriminant:
    //  < 0 has imaginary roots - no intersection
    //  = 0 has exact one root - we're at the exact edge of the sphere
    //  > 0 has two roots - the ray goes through the sphere
    const discriminant = b * b - 4 * a * c;

    if (discriminant < 0) {
        // Ray doesn't touch the sphere, return infinity.
        t1.* = std.math.inf(f32);
        t2.* = std.math.inf(f32);
        return;
    }

    t1.* = (-1 * b + @sqrt(discriminant)) / (2 * a);
    t2.* = (-1 * b - @sqrt(discriminant)) / (2 * a);
}

// Image memory! (width*height*3 bytes)
var img: [height][width]Color = .{.{Color{}} ** width} ** height;

// Takes coordinate system x,y points where 0,0 is the center of the image.
fn putPixel(x: i32, y: i32, c: Color) void {
    const img_x: u32 = @intCast(u32, (width / 2) + x);
    const img_y: u32 = @intCast(u32, (height / 2) - y);
    img[img_y][img_x] = c;
}

fn writeTga() !void {
    const file = try std.fs.cwd().createFile("foo.tga", .{ .read = true });
    defer file.close();
    // https://stackoverflow.com/a/49658800/695615
    // https://en.wikipedia.org/wiki/Truevision_TGA
    // http://paulbourke.net/dataformats/tga/
    // Note that all multi-byte values are little-endian.
    const header: [18]u8 = .{
        0, // 1  - No ID field (length of 0)
        0, // 2  - No color map
        2, // 3  - Uncompressed true-color image
        0, // 4  \
        0, // 5  |
        0, // 6  | Color map information (none)
        0, // 7  |
        0, // 8  /
        0, // 9  \ Origin X (16 bits)
        0, // 10 /
        0, // 11 \ Origin Y (16 bits)
        0, // 12 /
        width & 255, // 13 \ Width (px) mask last 8 bits
        (width >> 8) & 255, // 14 / Width (px) right shift and mask to get first 8 bits
        height & 255, // 15 \ Height (px) last
        (height >> 8) & 255, // 16 / Height (px) first
        24, // 17 - Bits per pixel (3 colors, 8 bits each)
        0b00100000, // 18 - Image descriptor (Bits 4,5 are origin. Set to "top left".)
    };

    try file.writeAll(&header);

    var w: u32 = 0; // Pixel counter per row (width)
    var h: u32 = 0; // Row counter (height)

    var out_buffer: [width * 3]u8 = undefined;

    while (h < height) : (h += 1) {
        w = 0;
        while (w < width) : (w += 1) {
            // IMPORTANT: TGA stores in BGR order, not RGB.
            out_buffer[w * 3] = img[h][w].b;
            out_buffer[w * 3 + 1] = img[h][w].g;
            out_buffer[w * 3 + 2] = img[h][w].r;
        }

        // Write buffered row
        try file.writeAll(&out_buffer);
    }
}

//
// Tests!
//
const expect = @import("std").testing.expect;

test "Vec3 minus" {
    const a = Vec3.new(1, 2, 3);
    const b = Vec3.new(3, 2, 1);
    const c = a.minus(b);
    try expect(c.x == -2.0);
    try expect(c.y == 0.0);
    try expect(c.z == 2.0);
}

test "Vec3 dot product" {
    const a = Vec3.new(1, 2, 3);
    const b = Vec3.new(3, 2, 1);
    const c = a.dot(b);
    try expect(c == 10);
}

test "Vec3 scale" {
    const a = Vec3.new(1, 2, 3);
    const c = a.scale(5);
    try expect(c.x == 5);
    try expect(c.y == 10);
    try expect(c.z == 15);
}
