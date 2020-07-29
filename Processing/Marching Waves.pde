import java.util.*;
import processing.svg.*;
import ch.bildspur.postfx.builder.*;
import ch.bildspur.postfx.pass.*;
import ch.bildspur.postfx.*;

PostFX fx;

String filename = "hand2.png";
PImage reference, tracemap;
PShape stencil;
PGraphics canvas;
int interval = 8;
int fieldNo = 0;
int currentFrame = 0;
int step;
int gridsize = 100;
color line = color(0);
color bg = color(255);
color from = #B77A64;
color to = #FFFF00;

float[][] SobelX = {{1, 0, -1}, {2, 0, -2}, {1, 0, -1}};
float[][] SobelY = {{-1, -2, -1}, {0, 0, 0}, {1, 2, 1}};
float[][] erode = {{1.0/9, 1.0/9, 1.0/9},
                   {1.0/9, 1.0/9, 1.0/9},
                   {1.0/9, 1.0/9, 1.0/9}};

PriorityQueue<Integer> narrowBand;
IntList origin, obstacles;
ArrayList<ArrayList<IntList>> snapshots;
ArrayList<IntList> paths;
PVector[] gradients, tangents;
boolean[] frozen, obstacleMap;
float[] T, F, D, tangentX, tangentY, curl, boundaryMap;

boolean render = false, pencil = false, record = false, run = true,
        showField = false, tSolved = false, fieldSolved = false, saved = false,
        overlay = false;

void setup(){
   fullScreen(P2D);
  //size(600, 600, P2D);
  canvasSetup();
  fx = new PostFX(this);
  // pixelDensity(2);
  // println(height);
}

void draw(){
  if(saved){
    delay(2000);
    saved = false;
  }
  fieldNo %= 3;
  currentFrame %= interval * 2;
  if(render){
    int s = narrowBand.size();
    for(int i = 0; i < s; i ++){
      if(!narrowBand.isEmpty()){
        updateMap();
      }
    }
    if(tSolved && !fieldSolved){
      for(int o: obstacles){
        for(int n: searchArea(o)){
          if(T[n] < T[o]) T[o] = T[n] + 1;
        }
      }
      for(int i = 0; i < pixels.length; i ++){
        float x = convolution(i, SobelX, T);
        float y = convolution(i, SobelY, T);
        PVector v = new PVector(x, y);
        gradients[i] = v;
        tangents[i] = new PVector(-y, x);
        v.rotate(HALF_PI);
        tangentX[i] = v.x;
        tangentY[i] = v.y;
      }
      boundaryMap = erode(boundaryMap);
      for(int i = 0; i < pixels.length; i ++) curl[i] = abs(convolution(i, SobelX, tangentY) - convolution(i, SobelY, tangentX));
      fieldSolved = true;
      // for(IntList i: snapshots.get(currentFrame)) trace(i.copy());
    }
    if(narrowBand.isEmpty()){
      tSolved = true;
      if(!showField){
        background(bg);
        loadPixels();
        animate(currentFrame);
        updatePixels();
        delay(50);
        if(run) currentFrame ++;
      }else{
        colorMode(HSB, 360, 255, 255);
        loadPixels();
        for(int i = 0; i < pixels.length; i ++){
          if(fieldNo == 0) pixels[i] = color(T[i] / 2);
          else if(fieldNo == 1) pixels[i] = color(degrees(gradients[i].heading() + PI), 255, 255);
          //else if(fieldNo == 2) pixels[i] = color(degrees(tangents[i].heading() + PI), 255, 255);
          else if(fieldNo == 2) pixels[i] = color(curl[i] * 10);
        }
        for(int i: obstacles) pixels[i] = color(0);
        // for(IntList i: paths){
        //   for(int j = 0; j < i.size(); j ++){
        //     pixels[i.get(j)] = color((j * 2) % 360, 255, 255);
        //   }
        // }
        // if(paths.isEmpty()){
        //   for(IntList i: snapshots.get(currentFrame)) trace(i.copy());
        // }
        updatePixels();
        // image(tracemap, 0, 0);
        colorMode(RGB, 255, 255, 255);
      }
    }
  }
  else{
    background(0, 0, 255);
    image(reference, 0, 0);
    canvas.beginDraw();
    canvas.stroke(255, 0, 0);
    canvas.strokeWeight(2);
    canvas.fill(255, 0, 0);
    if(keyPressed){
      if((key == 'r' || key == 'R')){
        for(int i = 0; i < 2; i ++) canvas.point(random(0, width), random(0, height));
      }
      else if((key == 'l' || key == 'L')){
        canvas.line(random(0, width), 0, random(0, width), height - 1);
        canvas.line(0, random(0, height), width - 1, random(0, height));
      }
      else if((key == 'g' || key == 'G')){
        int a = gridsize;
        //pushMatrix();
        //translate(width/2 - 150, -100);
        //rotate(radians(45));
        for(int i = a; i < height; i += a) canvas.line(0, i, width, i);
        for(int i = a; i < width; i += a) canvas.line(i, 0, i, height);
        //popMatrix();
      }
    }
    if(overlay){
      fill(255, 0, 0);
      shape(stencil, mouseX, mouseY);
      if(pencil) canvas.shape(stencil, mouseX, mouseY);
    }
    else if(pencil) canvas.line(pmouseX, pmouseY, mouseX, mouseY);
    canvas.endDraw();
    image(canvas, 0, 0);
  }

  //int bloom = 20;
  //if(fieldSolved){
  //  fx.render()
  //    //.sobel()
  //    .bloom(0.1, bloom, bloom)
  //    .compose();
  //}
}

void mousePressed(){pencil = true;}
void mouseReleased(){pencil = false;}

void keyReleased(){
  if(key == ENTER){
    render = !render;
    if(!render) canvasSetup();
    else{
      loadPixels();
      initialize();
    }
  }
  else if((key == 's' || key == 'S') && !run && fieldSolved) save(currentFrame);
  else if(key == 'f' || key == 'F') showField = !showField;
  else if(key == ' ' && !showField) run = !run;
  else if(key == 'z' || key == 'Z') canvas.clear();
  else if(key == 'k' || key == 'K') overlay = !overlay;
  else if(key == 'q' || key == 'Q') fieldNo ++;
  else if(key == CODED){
    if(keyCode == LEFT){
      if(currentFrame > 0) currentFrame --;
      else currentFrame = interval * 2 - 1;
    }
    else if(keyCode == RIGHT) currentFrame ++;
  }
}

color frameColor(int frame){
  float c = map(frame, 0, interval * 2, 0, 1);
  return lerpColor(from, to, c);
}

void animate(int frame){
  ArrayList<IntList> f = snapshots.get(frame);
  //color c;
  //if(colorIn) c = frameColor((frame + interval * 2 - (frameCount % (interval * 2))) % (interval * 2));
  //else c = line;
  for(IntList i: f){
    for(int j: i){
      pixels[j] = line;
    }
  }
}

void snapshot(float d){
  float t = (float)d;
  IntList i = new IntList();
  //background(bg);
  loadPixels();
  for(int c: narrowBand){
    //if(step % (interval * 2) == 0) pixels[c] = line;
    pixels[c] = line;
    i.append(c);
  }
  updatePixels();
  snapshots.get(step % (interval * 2)).add(i);
  step ++;
}

void canvasSetup(){
  canvas = createGraphics(width, height);
  tSolved = false;
  fieldSolved = false;
  colorMode(RGB, 255, 255, 255);
  step = 2;
  narrowBand = new PriorityQueue<Integer>(new cellComparator());
  origin = new IntList();
  obstacles = new IntList();
  snapshots = new ArrayList<ArrayList<IntList>>(interval * 2);
  for(int i = 0; i < interval * 2; i++){
    snapshots.add(new ArrayList<IntList>());
  }
  paths = new ArrayList<IntList>();
  frozen = new boolean[width * height];
  obstacleMap = new boolean[width * height];
  F = new float[width * height];
  T = new float[width * height];
  tangentX = new float[width * height];
  tangentY = new float[width * height];
  curl = new float[width * height];
  boundaryMap = new float[width * height];
  D = new float[width * height];
  gradients = new PVector[width * height];
  tangents = new PVector[width * height];
  reference = loadImage(filename);
  reference.resize(0, height);
  tracemap = createImage(width, height, RGB);
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

void initialize(){
  loadPixels();
  for(int i = 0; i < width * height; i ++){
    frozen[i] = false;
    obstacleMap[i] = false;
    color sample = reference.pixels[convertCoords(i, reference)];
    color sample2 = pixels[i];
    F[i] = map(brightness(sample), brightness(line), brightness(bg), 0.26, 2);
    if(blue(sample) > 230 && saturation(sample) > 230) obstacles.append(i);
    boundaryMap[i] = 1;
    if(red(sample2) > 230 && saturation(sample2) > 230){
      origin.append(i);
      frozen[i] = true;
      T[i] = 0;
    }
  }
  //colorMode(HSB);
  background(bg);
  for(int i = 0; i < width; i ++){
    obstacles.append(i);
    obstacles.append(width * height - i - 1);
  }
  for(int i = 0; i < height; i ++){
    obstacles.append(i * width);
    obstacles.append((i + 1) * width - 1);
  }
  for(int o: obstacles){
    frozen[o] = true;
    obstacleMap[o] = true;
    T[o] = Float.MAX_VALUE;
    boundaryMap[o] = 0;
  }
  for(int o: origin){
    for(int n: neighbors(o)){
      if(!frozen[n]){
        T[n] = 1/F[n];
        narrowBand.add(n);
      }
    }
  }
}

void updateMap(){
  int c = narrowBand.poll();
  frozen[c] = true;
  if(T[c] > float(step) / 2) snapshot(T[c]);
  for(int n: neighbors(c)){
      if(!frozen[n]){
        T[n] = solveEikonal(n);
        if(narrowBand.contains(n)) narrowBand.remove(n);
        narrowBand.add(n);
      }
    }
}

float solveEikonal(int u){
  float a, b, c, f = sq(F[u]);
  a = b = c = 0;
  for(int d = 0; d < 2; d ++){
    for(int j = 1; j > -2; j -= 2){
      float v = Float.MAX_VALUE;
      int n = getCell(d, u, j);
      if(n != -1 && frozen[n]){
        if(T[n] < v){
          v = (float)T[n];
        }
      }
      if(v < Float.MAX_VALUE  && !obstacleMap[n]){
        a += f;
        b -= 2 * f * v;
        c += sq(v) * f;
      }
    }
  }
  c -= 1;
  return solveQuadratic(a, b, c);
}

int getCell(int dimension, int c, int delta){
  int r = -1;
  if(dimension == 1){
    if((c % width) + 1 < width && (c % width) - 1 >= 0) r = c + delta;
  }
  else if(dimension == 0){
    if(c + width < T.length && c - width >= 0) r = c + (delta * width);
  }
  return r;
}

IntList neighbors(int c){
  IntList n = new IntList();
  if((c % width) + 1 < width) n.append(c + 1);
  if((c % width) - 1 >= 0) n.append(c - 1);
  if(c + width < T.length) n.append(c + width);
  if(c - width >= 0) n.append(c - width);
  return n;
}

float convolution(int sample, float[][] kernel, float[] image){
  int x = sample % width;
  int y = sample / width;
  float total = 0;
  int offset = 1;
  int matrixSize = 3;
  for (int i = 0; i < matrixSize; i++){
    for (int j= 0; j < matrixSize; j++){
      int xloc, yloc;
      if(x > 0 && x < width - 1) xloc = x - i + offset;
      else xloc = x;
      if(y > 0 && y < height - 1) yloc = y - j + offset;
      else yloc = y;
      total += image[loc(xloc, yloc)] * kernel[i][j];
    }
  }
  return total;
}

int loc(int x, int y){
  return x + y * width;
}

float solveQuadratic(float a, float b, float c){
  float det = sq(b) - (4 * a * c);
  return (-b + sqrt(det)) / (2 * a);
}

int convertCoords(int i, PImage to){
  int x = i % width;
  int y = i / width;
  return constrain(x + y * to.width, 0, to.width * to.height - 1);
}

void trace(IntList i){
  for(int k: i){
    if(boundaryMap[k] < 1.0) i.removeValue(k);
  }
  i.shuffle();
  while(i.size() > 0){
    boolean finishPath[] = {false, false};
    int[] seekers = new int[2];
    IntList[] sides = {new IntList(), new IntList()};
    seekers[0] = i.get(0);
    tracemap.pixels[seekers[0]] = color(255, 255, 0);
    i.removeValue(seekers[0]);
    sides[0].append(seekers[0]);
    seekers[1] = abs(closest(seekers[0], i, 1));
    tracemap.pixels[seekers[1]] = color(0, 255, 255);
    i.removeValue(seekers[1]);
    sides[1].append(seekers[1]);
    if(weightedDistance(seekers[1], closest(seekers[1], i, 1), 1) < 0.3) finishPath[1] = true;
    if(weightedDistance(seekers[0], closest(seekers[0], i, 0), 0) < 0.3) finishPath[0] = true;
    while(!(finishPath[0] && finishPath[1])){
      for(int j = 1; j > -1; j --){
        int c = closest(seekers[j], i, j);
        if(!finishPath[j] && c >= 0 && !sides[(j + 1) % 2].hasValue(c)){
          tracemap.pixels[c] = color(j * 255, 0, 255);
          sides[j].append(c);
          seekers[j] = c;
          i.removeValue(c);
        } else finishPath[j] = true;
      }
    }
    sides[1].reverse();
    IntList path = new IntList();
    path.append(sides[1]);
    path.append(sides[0]);
    if(path.size() > 5) paths.add(path);
   }
}

int closest(int a, IntList i, int dir){
  int c = -1;
  float dist = Float.MAX_VALUE;
  for(int b: searchArea(a)){
    if(i.hasValue(b) && weightedDistance(a, b, dir) < dist){
      c = b;
      dist = weightedDistance(a, b, dir);
    }
  }
  return c;
}

IntList searchArea(int c){
  IntList n = new IntList();
  for(int dx = -2; dx < 3; dx ++){
    for(int dy = -2; dy < 3; dy ++){
      int l = loc((c % width) + dx, (c / width) + dy);
      if(l >= 0 && l < width * height && l != c){
        n.append(l);
      }
    }
  }
  return n;
}

float weightedDistance(int a, int b, int dir){
  float d = solveEuclidean(a, b);
  PVector aV = toVector(a);
  PVector bV = toVector(b);
  PVector t = tangents[a].copy();
  t.rotate(HALF_PI);
  bV.sub(aV);
  float theta = PVector.angleBetween(bV, t);
  // theta = 1 - (abs(theta - HALF_PI) * 1.5 / PI);
  theta /= TWO_PI;
  if(dir == 0) theta += 0.5;
  else theta = 1 - theta;
  return d * theta;
}

float solveEuclidean(int a, int b){
  int x1 = a % width;
  int y1 = a / width;
  int x2 = b % width;
  int y2 = b / width;
  return sqrt(sq(x1 - x2) + sq(y1 - y2));
}

PVector toVector(int i){
  return new PVector(i % width, i / width);
}

IntList intersection(IntList a, IntList b){ //Largest list must go first
  IntList i = new IntList();
  for(int j: a){
    if(b.hasValue(j)) i.append(j);
  }
  return i;
}

public class cellComparator implements Comparator<Integer>{
  @Override
  public int compare(Integer a, Integer b){
    if(T[a] < T[b]) return -1;
    if(T[a] > T[b]) return 1;
    else return 0;
  }
}

IntList dither(IntList k){
  IntList path = k.copy();
  IntList d = new IntList();
  float c = 0;
  d.append(path.get(0));
  for(int j = 0; j < path.size(); j ++){
    int i = path.get(j);
    c += curl[i];
    if(c > 20){
      d.append(i);
      c = 0;
      j ++;
    }
  }
  return d;
}

float[] erode(float[] map){
  float[] n = new float[map.length];
  for(int i = 0; i < n.length; i ++) n[i] = floor(convolution(i, erode, map));
  return n;
}

void save(int frame){
  paths.clear();
  for(IntList i: snapshots.get(frame)) trace(i.copy());
  strokeWeight(1);
  noFill();
  String s = "Capture " + hour() + ":" + minute() + ":" + second() + ".svg";
  beginRecord(SVG, s);
  background(bg);
  for(int r = 0; r < paths.size(); r ++){
    IntList i = paths.get(r);
    IntList d = dither(i);
    beginShape();
    PVector v = toVector(d.get(0));
    vertex(v.x, v.y);
    for(int k = 1; k < d.size(); k ++){
      int a = d.get(k - 1);
      int b = d.get(k);
      float dist;
      dist = solveEuclidean(a, b) / 2;
      PVector tangentFront = gradients[b].copy();
      PVector tangentBack = gradients[a].copy();
      tangentFront.rotate(HALF_PI);
      tangentBack.rotate(HALF_PI);
      tangentFront.setMag(-dist);
      tangentBack.setMag(dist);
      PVector aloc = toVector(a);
      PVector bloc = toVector(b);
      tangentFront.add(bloc);
      tangentBack.add(aloc);
      bezierVertex(tangentBack.x, tangentBack.y, tangentFront.x, tangentFront.y, bloc.x, bloc.y);
    }
    stroke(line);
    noFill();
    PVector w = toVector(i.get(i.size() - 1));
    w.lerp(v, 0.5);
    if(boundaryMap[loc(round(w.x), round(w.y))] == 1.0 && PVector.dist(v, w) < 5) endShape(CLOSE);
    else endShape();
  }
  endRecord();
  println("recorded");
  saved = true;
}
