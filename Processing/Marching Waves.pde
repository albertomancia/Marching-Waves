import java.util.*;
import processing.svg.*;

String filename = "ellie.png";
PImage reference;
PShape stencil;
PGraphics canvas;
int interval = 14; //Try changing this!
int currentFrame = 0;
int step, count;
int gridsize = 100;
int res = 10;
color line = color(0);
color bg = color(255);

float epsilon = 0.1;
float threshold = 10;
float[][] erode = {{1.0/9, 1.0/9, 1.0/9},
                   {1.0/9, 1.0/9, 1.0/9},
                   {1.0/9, 1.0/9, 1.0/9}};

PriorityQueue<Integer> narrowBand;
IntList origin, obstacles;
ArrayList<ArrayList<IntList>> snapshots;
ArrayList<Segment> paths = new ArrayList<Segment>();
//Each pixel is referenced by a single integer, which can access each of these arrays:
boolean[] frozen, obstacleMap;
float[] T, F, boundaryMap;

boolean render = false, pencil = false, record = false, run = true,
  showField = false, tSolved = false, saved = false, overlay = false;

void setup() {
  fullScreen(P2D);
  //size(800, 600, P2D);
  canvasSetup();
  paths = new ArrayList<Segment>();
}

void draw(){
  if (saved) {
    delay(2000);
    saved = false;
  }
  currentFrame %= interval * 2;
  //Instructions for render mode
  if(render) {
    //For each step, process every pixel in the narrow band
    int s = narrowBand.size();
    for (int i = 0; i < s; i ++) {
      if (!narrowBand.isEmpty()) {
        updateMap();
      }
    }
    //Display and animate the lines
    if (narrowBand.isEmpty()) {
      if(!tSolved) boundaryMap = erode(boundaryMap);
      tSolved = true;
      if (!showField) {
        background(bg);
        loadPixels();
        animate(currentFrame);
        updatePixels();
        delay(25);
        if (run) currentFrame ++;
      } else {
        loadPixels();
        for (int i = 0; i < pixels.length; i ++) {
          pixels[i] = color(T[i] /4);
        }
        for (int i : obstacles) pixels[i] = color(0);
        updatePixels();
      }
    }
  }
  else drawingMode();
}

void mousePressed() {
  pencil = true;
}

void mouseReleased() {
  pencil = false;
}

//Key bindings
void keyReleased() {
  if (key == ENTER) {
    render = !render;
    if (!render) canvasSetup();
    else {
      loadPixels();
      initialize();
    }
  }
  else if ((key == 's' || key == 'S') && !run && tSolved) save();
  else if (key == 'f' || key == 'F') showField = !showField;
  else if (key == ' ' && !showField) run = !run;
  else if (key == 'z' || key == 'Z') canvas.clear();
  else if (key == 'k' || key == 'K') overlay = !overlay;
  else if (key == CODED) {
    if (keyCode == LEFT) {
      if (currentFrame > 0) currentFrame --;
      else currentFrame = interval * 2 - 1;
    } else if (keyCode == RIGHT) currentFrame ++;
  }
}

void animate(int frame) {
  ArrayList<IntList> f = snapshots.get(frame);
  for (IntList i : f) {
    for (int j : i) {
      pixels[j] = line;
    }
  }
}

void canvasSetup() {
  canvas = createGraphics(width, height);
  tSolved = false;
  colorMode(RGB, 255, 255, 255);
  step = 0;
  narrowBand = new PriorityQueue<Integer>(new cellComparator());
  //narrowBand = new LinkedList<Integer>();
  origin = new IntList();
  obstacles = new IntList();
  snapshots = new ArrayList<ArrayList<IntList>>(interval * 2);
  for (int i = 0; i < interval * 2; i++) {
    snapshots.add(new ArrayList<IntList>());
  }
  frozen = new boolean[width * height];
  obstacleMap = new boolean[width * height];
  F = new float[width * height];
  T = new float[width * height];
  boundaryMap = new float[width * height];
  reference = loadImage(filename);
  reference.resize(0, height);
  //Turn this on to create a noise pattern as the reference image
  //float increment = 0.02;
  //reference.loadPixels();
  //for(int i = 0; i < width * height; i ++){
  //  float n = noise((i % width) * increment, (i / width) * increment);
  //  n *= 255;
  //  reference.pixels[i] = color(2 * (n - 128) + 128);
  //}
  //reference.updatePixels();
  //background(100);
  background(0, 0, 255);
  image(reference, 0, 0);
}

public class cellComparator implements Comparator<Integer> {
  @Override
    public int compare(Integer a, Integer b) {
    if (T[a] < T[b]) return -1;
    if (T[a] > T[b]) return 1;
    else return 0;
  }
}

void initialize() {
  loadPixels();
  for (int i = 0; i < width * height; i ++) {
    T[i] = Float.MAX_VALUE;
    frozen[i] = false;
    obstacleMap[i] = false;
    color sample = reference.pixels[convertCoords(i, reference)];
    color sample2 = pixels[i];
    F[i] = map(brightness(sample), brightness(line), brightness(bg), 0.13, 1);
    if (blue(sample) > 230 && saturation(sample) > 230) obstacles.append(i);
    boundaryMap[i] = 1;
    if (red(sample2) > 230 && saturation(sample2) > 230) {
      origin.append(i);
      frozen[i] = true;
      T[i] = 0;
    }
  }
  background(bg);
  for (int i = 0; i < width; i ++) {
    obstacles.append(i);
    obstacles.append(width * height - i - 1);
  }
  for (int i = 0; i < height; i ++) {
    obstacles.append(i * width);
    obstacles.append((i + 1) * width - 1);
  }
  for (int o : obstacles) {
    frozen[o] = true;
    obstacleMap[o] = true;
    T[o] = Float.MAX_VALUE;
    boundaryMap[o] = 0;
  }
  for (int o : origin) {
    for (int n : neighbors(o)) {
      if (!frozen[n] && !narrowBand.contains(n)) {
        T[n] = 1/F[n];
        narrowBand.add(n);
      }
    }
  }
}

int convertCoords(int i, PImage to) {
  int x = i % width;
  int y = i / width;
  return constrain(x + y * to.width, 0, to.width * to.height - 1);
}

void updateMap() {
  int c = narrowBand.poll();
  frozen[c] = true;
  if (T[c] > float(step) / 2) snapshot();
  for (int n : neighbors(c)) {
    if (!frozen[n]) {
      T[n] = solveEikonal(n);
      if (narrowBand.contains(n)) narrowBand.remove(n);
      narrowBand.add(n);
    }
  }
}

IntList neighbors(int c) {
  IntList n = new IntList();
  if ((c % width) + 1 < width) n.append(c + 1);
  if ((c % width) - 1 >= 0) n.append(c - 1);
  if (c + width < T.length) n.append(c + width);
  if (c - width >= 0) n.append(c - width);
  return n;
}

void snapshot() {
  IntList i = new IntList();
  loadPixels();
  for (int c : narrowBand) {
    if(step % (interval * 2) == 0) pixels[c] = line;
    i.append(c);
  }
  updatePixels();
  snapshots.get(step % (interval * 2)).add(i);
  step ++;
}

float solveEikonal(int u) {
  float a, b, c, f = sq(F[u]);
  a = b = c = 0;
  for (int d = 0; d < 2; d ++) {
    float v = Float.MAX_VALUE;
    for (int j = 1; j > -2; j -= 2) {
      int n = getCell(d, u, j);
      if (n != -1 && frozen[n] && T[n] < v && !obstacleMap[n]) {
        v = (float)T[n];
      }
    }
    if (v < Float.MAX_VALUE) {
      a += f;
      b -= 2 * f * v;
      c += sq(v) * f;
    }
  }
  c -= 1;
  return solveQuadratic(a, b, c);
}

int getCell(int dimension, int c, int delta) {
  int r = -1;
  if (dimension == 1) {
    if ((c % width) + 1 < width && (c % width) - 1 >= 0) r = c + delta;
  } else if (dimension == 0) {
    if (c + width < T.length && c - width >= 0) r = c + (delta * width);
  }
  return r;
}

float solveQuadratic(float a, float b, float c) {
  float det = sq(b) - (4 * a * c);
  return (-b + sqrt(det)) / (2 * a);
}

float[] erode(float[] map){
  float[] n = new float[map.length];
  for(int i = 0; i < n.length; i ++) n[i] = floor(convolution(i, erode, map));
  return n;
}

float convolution(int sample, float[][] kernel, float[] image) {
  int x = sample % width;
  int y = sample / width;
  float total = 0;
  int offset = 1;
  int matrixSize = 3;
  for (int i = 0; i < matrixSize; i++) {
    for (int j= 0; j < matrixSize; j++) {
      int xloc, yloc;
      if (x > 0 && x < width - 1) xloc = x - i + offset;
      else xloc = x;
      if (y > 0 && y < height - 1) yloc = y - j + offset;
      else yloc = y;
      total += image[loc(xloc, yloc)] * kernel[i][j];
    }
  }
  return total;
}

int loc(int x, int y) {
  return x + y * width;
}

void drawingMode(){
  background(0, 0, 255);
  image(reference, 0, 0);
  canvas.beginDraw();
  canvas.stroke(255, 0, 0);
  canvas.strokeWeight(2);
  canvas.fill(255, 0, 0);
  //Automatic drawing functions
  if (keyPressed) {
    if ((key == 'r' || key == 'R')) {
      for (int i = 0; i < 2; i ++) canvas.point(random(0, width), random(0, height));
    } else if ((key == 'l' || key == 'L')) {
      canvas.line(random(0, width), 0, random(0, width), height - 1);
      canvas.line(0, random(0, height), width - 1, random(0, height));
    } else if ((key == 'g' || key == 'G')) {
      int a = gridsize;
      //Turn this on to rotate the grid
      //pushMatrix();
      //translate(width/2 - 150, -100);
      //rotate(radians(45));
      for (int i = a; i < height; i += a) canvas.line(0, i, width, i);
      for (int i = a; i < width; i += a) canvas.line(i, 0, i, height);
      //popMatrix();
    }
  } else if (pencil) canvas.line(pmouseX, pmouseY, mouseX, mouseY);
  canvas.endDraw();
  image(canvas, 0, 0);
}

//SVG Export

void createLines(){
  noFill();
  strokeWeight(2);
  stroke(255, 0, 0);
  for (int i = 0; i < width-1; i ++) {
    for (int j = 0; j < height-1; j ++) {
      float x = i;
      float y = j;
      float[] v = {T[loc(i, j)], T[loc(i+1, j)], T[loc(i+1, j+1)], T[loc(i, j+1)]};
      Point a = new Point(x + lerp(v[0], v[1]), y                   );
      Point b = new Point(x + 1,                y + lerp(v[1], v[2]));
      Point c = new Point(x + lerp(v[3], v[2]), y + 1               );
      Point d = new Point(x,                    y + lerp(v[0], v[3]));
      int state = getState(signBit(v[0] - threshold), signBit(v[1] - threshold),
                           signBit(v[2] - threshold), signBit(v[3] - threshold));
      switch (state){
      case 1:
        addSegment(c, d);
        break;
      case 2:
        addSegment(b, c);
        break;
      case 3:
        addSegment(b, d);
        break;
      case 4:
        addSegment(a, b);
        break;
      case 5:
        addSegment(a, b);
        addSegment(c, d);
        break;
      case 6:
        addSegment(a, c);
        break;
      case 7:
        addSegment(a, d);
        break;
      case 8:
        addSegment(d, a);
        break;
      case 9:
        addSegment(c, a);
        break;
      case 10:
        addSegment(d, a);
        addSegment(b, c);
        break;
      case 11:
        addSegment(b, a);
        break;
      case 12:
        addSegment(d, b);
        break;
      case 13:
        addSegment(c, b);
        break;
      case 14:
        addSegment(d, c);
        break;
      }
    }
  }
  count = paths.size();
}

class Segment{
  Point[] points;
  Segment(Point a, Point b){
    Point[] h = {a, b};
    points = h;
  }
  Segment(Segment a, Segment b){
    Point[] h = (Point[])concat(shorten(a.points), b.points);
    points = h;
  }
}

class Point{
  float x, y;
  Point(float xpos, float ypos){
    x = xpos;
    y = ypos;
  }
}

int getState(int a, int b, int c, int d) {
  return a * 8 + b * 4  + c * 2 + d * 1;
}

int signBit(float f){
  return (Float.floatToIntBits(f) >> 31) + 1;
}

void addSegment(Point start, Point end) {
  if(boundaryMap[loc(round(start.x), round(start.y))] == 1){
    paths.add(new Segment(start, end));
  }
}

float lerp(float a, float b){
    return (a - threshold) / (a - b);
}

void sortstep(ArrayList<Segment> list){
 count --;
 for(int i = 0; i < list.size() - 1; i++){
   if(match(list.get(i), list.get(i+1))){
     combine(list, i);
     count = list.size();
   }
   if(i >= list.size() - 1) break;
   if(!match(list.get(i), list.get(i+1))){
     swap(list, i);
   }
 }
}

boolean match(Segment a, Segment b){
  int n = a.points.length - 1;
  return equal(a.points[n], b.points[0]);
}

boolean equal(Point v1, Point v2){
  return v1.x == v2.x && v1.y == v2.y;
}

void swap(ArrayList<Segment> list, int index){
  Segment a = list.get(index);
  Segment b = list.get(index + 1);
  list.set(index, b);
  list.set(index + 1, a);
}

void combine(ArrayList<Segment> list, int index){
  Segment a = list.get(index);
  Segment b = list.get(index + 1);
  Segment combinedSegment = new Segment(a, b);
  list.remove(a);
  list.remove(b);
  list.add(index, combinedSegment);
}

void render(ArrayList<Segment> list){
  noFill();
  strokeWeight(1);
  stroke(line);
  for(Segment s: list){
    beginShape();
    for(Point p: s.points){
      vertex(p.x, p.y);
    }
    endShape();
  }
}

void save(){
  String s = "Capture " + hour() + ":" + minute() + ":" + second() + ".svg";
  beginRecord(SVG, s);
  background(bg);
  for(int i = currentFrame / 2; i < step; i += interval){
    threshold = i;
    createLines();
    while(count > 1){
      sortstep(paths);
    }
    render(paths);
    paths.clear();
  }
  endRecord();
  println("recorded");
  saved = true;
}
