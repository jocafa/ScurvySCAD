// Lines ----------------------------------------------------------------------
function linearBasis(a, b, t) =
  (1-t)*a +
  t*b;

function linearPoint(a, b, t) = [
  linearBasis(a[0], b[0], t),
  linearBasis(a[1], b[1], t)
];

// Quadratic Bezier Curves ----------------------------------------------------
function quadraticBezierBasis(a, b, c, t) =
  pow(1-t, 2)*a +
  2*(1-t)*t*b +
  pow(t,2)*c;

function quadraticBezierPoint(a, b, c, t) = [
  quadraticBezierBasis(a[0], b[0], c[0], t),
  quadraticBezierBasis(a[1], b[1], c[1], t)
];

module quadraticBezierShape(a, b, c, n=10, convex=true) {
  m = convex ? [(a[0] + c[0])/2, (a[1] + c[1])/2] : b;
  union() for (i = [0 : n-1]) {
    polygon(
      points = [
        m,
        quadraticBezierPoint(a, b, c, (i+1)/n),
        quadraticBezierPoint(a, b, c, i/n),
      ]
    );
  }
}

// Cubic Bezier Curves --------------------------------------------------------
function cubicBezierBasis(a, b, c, d, t) =
  pow(1-t,3)*a +
  3*pow(1-t,2)*t*b +
  3*(1-t)*pow(t,2)*c +
  pow(t,3)*d;

function cubicBezierPoint(a, b, c, d, t) = [
  cubicBezierBasis(a[0], b[0], c[0], d[0], t),
  cubicBezierBasis(a[1], b[1], c[1], d[1], t)
];

// Elliptical Arcs ------------------------------------------------------------
function ellipticalAngle(u, v) =
  sign(u[0]*v[1] - u[1]*v[0]) *
  acos(u[0]*v[0] + u[1]*v[1]) / (
    sqrt(u[0]*u[0] + u[1]*u[1]) *
    sqrt(v[0]*v[0] + v[1]*v[1])
  );

function ellipticalArcPoint(center, radii, rotation, angle) = [
  [cos(rotation), -sin(rotation)],
  [sin(rotation), cos(rotation)]
] * [
  radii[0] * cos(angle),
  radii[1] * sin(angle)
] + [
  center[0],
  center[1]
];

module ellipticalArcShape(center, radii, rotation, theta1, deltaTheta, n=10) {
  union() for (i = [0 : n-1]) {
    polygon(
      points = [
        center,
        ellipticalArcPoint(center, radii, rotation, theta1 + (deltaTheta * (i+1) / n)),
        ellipticalArcPoint(center, radii, rotation, theta1 + (deltaTheta * i / n))
      ]
    );
  }
}

module ellipticalArcShapeSVG(start, radii, rotation=0, large, sweep, end, n=10) {
  x1 = start[0];
  y1 = start[1];
  x2 = end[0];
  y2 = end[1];

  midpoint = [(x1 + x2) / 2, (y1 + y2) / 2];

  phi = rotation % 360;

  startPrime = [
    [cos(phi), sin(phi)],
    [-sin(phi), cos(phi)]
  ] * [
    (x1 - x2) / 2,
    (y1 - y2) / 2
  ];

  x1p = startPrime[0];
  y1p = startPrime[1];

  rxp = abs(radii[0]);
  ryp = abs(radii[1]);

  lamda = (x1p*x1p/rxp*rxp) + (y1p*y1p/ryp*ryp);

  rx = lamda <= 1 ? rxp : sqrt(lamda)*rxp;
  ry = lamda <= 1 ? ryp : sqrt(lamda)*ryp;

  rx = rxp;
  ry = ryp;

  centerPrime =
    (large == sweep ? -1 : 1) *
    sqrt(
      (rx*rx * ry*ry - rx*rx * y1p*y1p - ry*ry * x1p*x1p)
      /
      (rx*rx * y1p*y1p + ry*ry * x1p*x1p)
    ) * [
      rx * y1p / ry,
      -(ry * x1p / rx)
    ];

  cxp = centerPrime[0];
  cyp = centerPrime[1];

  center = [
    [cos(phi), -sin(phi)],
    [sin(phi), cos(phi)]
  ] * [
    cxp,
    cyp
  ] + [
    (x1 + x2) / 2,
    (y1 + y2) / 2
  ];

  theta1= ellipticalAngle(
    [1, 0],
    [(x1p - cxp) / rx, (y1p - cyp) / ry]
  );

  deltaThetaPrime = ellipticalAngle(
    [(x1p - cxp) / rx, (y1p - cyp) / ry],
    [(-x1p - cxp) / rx, (-y1p - cyp) / ry]
  ) % 360;

  //theta1 = theta1prime;

  deltaTheta =
    (sweep == 0 && deltaThetaPrime > 0)
    ? deltaThetaPrime - 360
    : (sweep == 1 && deltaThetaPrime < 0)
      ? deltaThetaPrime + 360
      : deltaThetaPrime;

  //echo("center", center, "theta1", theta1, "sweep", sweep, "deltaTheta", deltaTheta, "rx", rx, "ry", ry);

  union() for (i = [0 : n-1]) {
    polygon(
      points = [
        midpoint,
        ellipticalArcPoint(center, [rx, ry], rotation, theta1 + (deltaTheta * (i+1) / n)),
        ellipticalArcPoint(center, [rx, ry], rotation, theta1 + (deltaTheta * i / n))
      ]
    );
  }

  //translate(center) cylinder(r = 0.5, h = 1);
}

color("Purple") {
  linear_extrude(height = 1) quadraticBezierShape([0, 0], [5, 10], [10, 0], n=20);
  translate([10, 0]) linear_extrude(height = 1) quadraticBezierShape([0, 0], [5, 10], [10, 0], n=20, convex=false);
}

startPoint = [0, 5];
endPoint = [-5, 0];
radii = [4, 8];
rotation = 0;

color("Black") translate(startPoint) cube(0.25, center=true);
color("Black") translate(endPoint) cube(0.25, center=true);

//color("Orange") translate([0, 0, -10]) ellipticalArcShape([0, 0], [20, 10], 80, 0, 90, n=10);

color("Red", 0.2)    translate([0, 0, 0]) linear_extrude(height = 1) ellipticalArcShapeSVG(startPoint, radii, rotation, 0, 0, endPoint, n=10);
color("Green", 0.2)  translate([0, 0, 2]) linear_extrude(height = 1) ellipticalArcShapeSVG(startPoint, radii, rotation, 0, 1, endPoint, n=10);
color("Blue", 0.2)   translate([0, 0, 4]) linear_extrude(height = 1) ellipticalArcShapeSVG(startPoint, radii, rotation, 1, 0, endPoint, n=10);
color("Yellow", 0.2) translate([0, 0, 6]) linear_extrude(height = 1) ellipticalArcShapeSVG(startPoint, radii, rotation, 1, 1, endPoint, n=10);

