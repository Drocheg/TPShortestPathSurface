import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.PriorityQueue;

import Jama.Matrix;
import Jcg.geometry.Point_3;
import Jcg.geometry.Vector_3;
import Jcg.polyhedron.Halfedge;
import Jcg.polyhedron.Polyhedron_3;
import Jcg.polyhedron.Vertex;

public class ExactGeodesics extends ShortestPathCalculator{

	static final double EPSILON = 0.0001;
	public ExactGeodesics(Polyhedron_3<Point_3> polyhedron3d) {
		super(polyhedron3d);
	}

	@Override
	public ArrayList<Point_3> calculatesShortestPath(Vertex<Point_3> source, Vertex<Point_3> destination) {
		// Step1: Use Dijkstra to compute upper-bound distance Ust(Dijkstra)
		// TODO test it for open forms? Or not?
		double ustDijkstra = bidirectionalDijkstra(source, destination);
		System.out.println("dijkstra: " + ustDijkstra);
		
		System.out.println("Source Vertex: " + source);
		System.out.println("Destination Vertex: " + destination);
		
		// Step2: Use approximation search to compute lowest bound distance function Lt

		// TODO
		// We can use for now as lowerBound function
		double lowestBoundDistance = lowestBoundDistance(source.getPoint(), destination.getPoint());
		System.out.println("lowestBoundDistance between source and destination: "+ lowestBoundDistance);
		
		
		
		// Step3: Use backtracking the create a tighter upper-bound Ust(backtrace) 
		// TODO. For now we still use Dijkstra
		
		
		// Step4: Do exact search to compute exact distance
		// TODO find why result is wrong (sometimes it is outside of bounds, sometimes even < lowerBound (imposible))
		// posible error for upperbound is that we end when arriving to first edge but other paths could exist.
		// Maybe continue until it is imposible, ingoring windows with lowest bound > current result.
		double upperBound = ustDijkstra;
		Map<Halfedge<Point_3>, ArrayList<Window>> mapWindowsEdge = new HashMap<Halfedge<Point_3>, ArrayList<Window>>();
		double shortestDistance = calculateShortestDistance(source, destination, upperBound, mapWindowsEdge);
		System.out.println("lowestBoundDistance between source and destination: "+ lowestBoundDistance);
		System.out.println("shortestDistance from vertex: " + source + "to vertex: " + destination + " is: " + shortestDistance);
		System.out.println("Dijkstra: " + ustDijkstra);
		System.out.println(lowestBoundDistance + " <= " + shortestDistance + " <= " + ustDijkstra);
		if(lowestBoundDistance>shortestDistance || shortestDistance > ustDijkstra) {
			//throw new IllegalStateException("It is not respecting the bounds");
			System.out.println("IT IS NOT RESPECTING BOUNDS!");
		}
		// Step5: Apply backtracking to find the shortest path.
		// TODO
		return backTracking(mapWindowsEdge, source, destination, shortestDistance);
	}
	
	
	private ArrayList<Point_3> backTracking(Map<Halfedge<Point_3>, ArrayList<Window>> mapWindowsEdge, Vertex<Point_3> source,
			Vertex<Point_3> destination, double shortestDistance) {
		// TODO Auto-generated method stub
		Map<Point_3, Halfedge<Point_3>> mapPointEdge = new HashMap();
		ArrayList<Point_3> listOfPoints = new ArrayList();
		Halfedge<Point_3> inicialHalfEdge = destination.getHalfedge(); // edge that point to the vertex
		Window bestWindow = findBestWindow(inicialHalfEdge, shortestDistance, mapWindowsEdge);
		
		
		Point_3 currentPoint = destination.getPoint();
		listOfPoints.add(currentPoint);
		Window currentWindow = bestWindow;
		mapPointEdge.put(currentPoint, currentWindow.edge);
		int iteration = 0;
		while(true) {
			System.out.println("Backpropagation iteration: " + iteration++);
			Point_3 p1_3D = currentWindow.edge.vertex.getPoint();
			Point_3 p0_3D = currentWindow.edge.opposite.vertex.getPoint();
			Point_3 p2_3D = currentWindow.edge.next.vertex.getPoint();
			double distanceToEdgeOrigin = euclidianDistance(currentPoint,  p0_3D);
			Point_3 p = new Point_3(distanceToEdgeOrigin,0,0);
			Point_3 vS = currentWindow.calculateSourceInXWithP0in000();
			if(isAlmostSamePoint(p,vS)){
				if(isAlmostSamePoint(currentPoint, source.getPoint())){
					break; //We are on the source
				}else {
					//We are on a source that is not destination
					if(!mapPointEdge.containsKey(currentPoint)) {
						throw new IllegalStateException("there is no edge with this point in the map.");
					}
					inicialHalfEdge = mapPointEdge.get(currentPoint); // edge that point to the vertex
					if(!isAlmostSamePoint(currentPoint, inicialHalfEdge.getVertex().getPoint())) {
						inicialHalfEdge = inicialHalfEdge.opposite;
						if(!isAlmostSamePoint(currentPoint, inicialHalfEdge.getVertex().getPoint())) throw new IllegalArgumentException("The edge do not contain the new source");
					}
					bestWindow = findBestWindow(inicialHalfEdge, currentWindow.sigma, mapWindowsEdge);
					currentWindow = bestWindow;
					//Do not change point
				}
			}else {
				//we are in a normal case.
				//calculate line between point and source and also lines of other 2 segments.
				Point_3 p2 = calculatePointWithDistances(size(currentWindow.edge), size(currentWindow.edge.prev), size(currentWindow.edge.next));
				Point_3 p0 = new Point_3(0,0,0);
				Point_3 p1 = new Point_3(size(currentWindow.edge),0,0);
				
				Vector_3 vp1_p2 = (Vector_3) p1.minus(p2);
				Vector_3 vp_s = (Vector_3) p.minus(vS);
				Vector_3 vp1_p = (Vector_3) p1.minus(p);
				
				Vector_3 vp2_p0 = (Vector_3) p2.minus(p0);
				Vector_3 vp2_p = (Vector_3) p2.minus(p);

				Point_3[] points = new Point_3[] {p1,p0,p2,p,vS,vS};
				
				double[] t = new double[2];
				if(isAlmostZero(vS.y)) {
					if(isAlmostZero(vS.x-p.x)) throw new IllegalStateException("already hanle this case");
					if(vS.x>p.x) {
						t[0]=0;
						t[1]=0;
					}else {
						t[0]=1;
						t[1]=1;
					}
				}else {
					calculateT(t, vp1_p, vp2_p, vp_s, vp1_p2, vp2_p0,points, 0);
				}
				Point_3 newPoint = null;
				if(isAlmostZero(t[0]-1) && isAlmostZero(t[1])){
					//it cut on the vertex. Might have problems
					newPoint = p2_3D;
					if(!isAlmostSamePoint(p2, vS)) {
						System.out.println("It pass by the vertex but it is not a source");
						// Test firstLine
						Halfedge<Point_3> newEdge = currentWindow.edge.next.opposite;
						double distanceToStartOfEdge = size(newEdge)*(1-t[0]);
						Window mightBeWindow = null;
						Window firstPosibleWindow = findPosibleWindow(mapWindowsEdge, newEdge, distanceToStartOfEdge);
						// Test secondLine
						newEdge = currentWindow.edge.prev.opposite;
						distanceToStartOfEdge = size(newEdge)*(1-t[1]);
						Window secondPosibleWindow = findPosibleWindow(mapWindowsEdge, newEdge, distanceToStartOfEdge);
						// CopyfromBellow
						if(firstPosibleWindow == null && secondPosibleWindow == null) throw new IllegalStateException("None window found with point. Problem backtracking");
						if(firstPosibleWindow == null) currentWindow = secondPosibleWindow;
						else if(secondPosibleWindow == null) currentWindow = firstPosibleWindow;
						else currentWindow = firstPosibleWindow.d0<secondPosibleWindow.d1?firstPosibleWindow:secondPosibleWindow;
					}else {
						//change window to a fake one for the algorithm
						Halfedge<Point_3> edgeForFakeWindow = currentWindow.edge.next;
						currentWindow = new Window(0, size(edgeForFakeWindow), size(edgeForFakeWindow), 0, currentWindow.sigma, false, edgeForFakeWindow, 0);
					}
				}else {
					Halfedge<Point_3> newEdge = null;
					double distanceToStartOfEdge = 0;
					if(t[0]<=1-EPSILON) {
						//it cut the first line
						newPoint = Point_3.linearCombination(new Point_3[]{p1_3D,p2_3D},new Double[] {1-t[0],t[0]});
						newEdge = currentWindow.edge.next.opposite;
						distanceToStartOfEdge = size(newEdge)*(1-t[0]);
					}else {
						if(t[1]>=EPSILON) {
							//it cut the second line
							newPoint = Point_3.linearCombination(new Point_3[]{p2_3D,p0_3D},new Double[] {1-t[1],t[1]});
							newEdge = currentWindow.edge.prev.opposite;
							distanceToStartOfEdge = size(newEdge)*(1-t[1]);
						}else {
							throw new IllegalStateException("bad calculation error for the cut in backpropagation");
						}
					}
					//change to new window.
					if(!mapWindowsEdge.containsKey(newEdge)) {
						throw new IllegalStateException("there was a problem calculating the backtracking. Non existing prev window");
					}
					ArrayList<Window> posibleNewWindows = mapWindowsEdge.get(newEdge);
					boolean foundWindow = false;
					for(Window posibleW: posibleNewWindows) {
						if(posibleW.b0-EPSILON < distanceToStartOfEdge && posibleW.b1+EPSILON>distanceToStartOfEdge) {
							if(foundWindow) {
								throw new IllegalStateException("Can't find two windows in same edge with same point");
							}
							currentWindow = posibleW;
							foundWindow = true;
						}
					}
					if(!foundWindow) {
						throw new IllegalStateException("None window found with point. Problem backtracking");
					}
				}
				currentPoint = newPoint;
				listOfPoints.add(currentPoint);
				mapPointEdge.put(currentPoint, currentWindow.edge);
			}
		}
		return listOfPoints;
	}


	private Window findPosibleWindow(Map<Halfedge<Point_3>, ArrayList<Window>> mapWindowsEdge,
			Halfedge<Point_3> newEdge, double distanceToStartOfEdge) {
		Window result = null;
		if(mapWindowsEdge.containsKey(newEdge)) {
			boolean foundWindow = false;
			ArrayList<Window> posibleNewWindows = mapWindowsEdge.get(newEdge);
			for(Window posibleW: posibleNewWindows) {
				if(posibleW.b0-EPSILON < distanceToStartOfEdge && posibleW.b1+EPSILON>distanceToStartOfEdge) {
					if(foundWindow) {
						throw new IllegalStateException("Can't find two windows in same edge with same point");
					}
					result = posibleW;
					foundWindow = true;
				}
			}
		}
		return result;
	}

	private Window findBestWindow(Halfedge<Point_3> inicialHalfEdge, double currentShortestDistance, Map<Halfedge<Point_3>, ArrayList<Window>> mapWindowsEdge) {
		Window bestWindow = null;
		Halfedge<Point_3> currentHalfEdge = inicialHalfEdge;	
		boolean stateNext = true;
		do {
			if(mapWindowsEdge.containsKey(currentHalfEdge)) {
				ArrayList<Window> windowArrray = mapWindowsEdge.get(currentHalfEdge);
				for(Window w: windowArrray) {
					if(!isAlmostZero(currentShortestDistance-w.sigma)) {
						if((isAlmostZero(w.d1+w.sigma-currentShortestDistance) && isAlmostZero(w.b1-size(currentHalfEdge)) && stateNext) || (isAlmostZero(w.d0+w.sigma-currentShortestDistance)&& isAlmostZero(w.b0)) && !stateNext) { 
							bestWindow = w;
						}
					}else {
						System.out.println("sigma==shortestDistance");
					}
				}
			}
			if(stateNext) currentHalfEdge = currentHalfEdge.next;
			else currentHalfEdge = currentHalfEdge.opposite;
			stateNext = !stateNext;
		}while(inicialHalfEdge!=currentHalfEdge);
		if(bestWindow==null) throw new IllegalArgumentException("there was no edge that touched the destination");
		return bestWindow;
	}

	public Point_3 calculatePointWithDistances(double d, double d0, double d1) {
			if(isAlmostZero(d-d0-d1)) {
				return new Point_3(d0,0,0);
			}
			if(isAlmostZero(d0-d-d1)) {
				return new Point_3(d0,0,0);
			}
			if(isAlmostZero(d1-d-d0)) {
				return new Point_3(-d0,0,0);
			}
			if(isAlmostZero(d)) throw new IllegalStateException("can't calculate point if sizeEdge==0");
			double a = (d0*d0-d1*d1+d*d)/(2*d);
			double h = Math.sqrt(d0*d0-a*a);
			if(Double.isNaN(h)) {
				throw new IllegalStateException("h is NaN in calculating point");
			}
			Point_3 pResult = new Point_3(a,h,0);
			return pResult;
		}

	private boolean isAlmostSamePoint(Point_3 p1, Point_3 p2) {
		return isAlmostZero(p1.distanceFrom(p2).doubleValue());
	}


	//Exact algorithm
	public static class Window implements Comparable<Window>{
		double b0, b1, d0, d1, sigma, pqCriteria;  //b0,b1 distance along the edge. d0,d1 distance to source. Sigma geodesic distance of source to Vs
		boolean tau; // side of the edge on which the source is
		Halfedge<Point_3> edge;
		public Window(double b0, double b1, double d0, double d1, double sigma, boolean tau, Halfedge<Point_3> edge, double pqCriteria) {
			super();
			this.b0 = b0;
			this.b1 = b1;
			this.d0 = d0;
			this.d1 = d1;
			this.sigma = sigma;
			this.tau = tau;
			this.edge = edge;
			this.pqCriteria = pqCriteria;
			if(Double.isNaN(pqCriteria)) {
				throw new IllegalArgumentException("pqCriteria can't be NaN");
			}
			Point_3 source = calculateSourceInXWithP0in000();
			if(!isAlmostZero(euclidianDistance(source, new Point_3(b0,0,0))-d0)) {
				throw new IllegalArgumentException("argument are not concording");
			}
			if(!isAlmostZero(euclidianDistance(source, new Point_3(b1,0,0))-d1)) {
				throw new IllegalArgumentException("argument are not concording");
			}
			System.out.println("New: " + toString());
			// TODO delete System.out.println("PQCriteria:  " + pqCriteria + " sigma: " + sigma);
		}
		
		
		@Override
		public String toString() {
			double pre = 100.0;
			return "Window [edge=" + edge + ", b0=" + Math.round(b0*pre)/pre + ", b1=" +  Math.round(b1*pre)/pre + ", d0=" +  Math.round(d0*pre)/pre + ", d1=" +  Math.round(d1*pre)/pre + ", sigma=" +  Math.round(sigma*pre)/pre
					+ ", pqCriteria=" + pqCriteria + ", tau=" + tau + "]";
		}


		@Override
		public int compareTo(Window w) {
			if(w == null) throw new IllegalArgumentException();
			double diff = pqCriteria-w.pqCriteria;
			if(diff>0) return 1;
			if(diff<0) return -1;
			return 0;
		}
		public ArrayList<Window> propagate(Vertex<Point_3> source3D) {
			ArrayList<Window> windows = new ArrayList<Window>();
			if(d0==0) {
				Halfedge<Point_3> newEdge = edge.opposite.next;
				Window w = new Window(0, size(newEdge), 0, size(newEdge), sigma, tau, newEdge, sigma);
				windows.add(w);
				Halfedge<Point_3> newEdge2  = newEdge.next;
				Window w2 = new Window(0, size(newEdge2), size(newEdge), size(edge), sigma, tau, newEdge2, sigma+Math.min(size(newEdge), size(edge)));
				windows.add(w2);
				System.out.println("d0==0: " + "first w: " + w + "second w: " + w2);
				return windows;
			}
			if(d1==0) {
				Halfedge<Point_3> newEdge = edge.opposite.prev;
				Window w = new Window(0, size(newEdge),  size(newEdge), 0, sigma, tau, newEdge, sigma);
				windows.add(w);
				Halfedge<Point_3> newEdge2  = newEdge.prev;
				Window w2 = new Window(0, size(newEdge2), size(edge), size(newEdge), sigma, tau, newEdge2, sigma+Math.min(size(newEdge), size(edge)));
				windows.add(w2);
				System.out.println("d1==0: " + "first w: " + w + "second w: " + w2);
				return windows;
			}
			
			if(tau == true) throw new IllegalStateException(); //tau is always false because of how the windows are propagated. To demonstrate.
			Point_3 p1 = new Point_3(edge.vertex.getPoint()); //end of edge
			Point_3 p0 = new Point_3(edge.opposite.vertex.getPoint()); //start of edge
			Point_3 p2 = new Point_3(edge.next.vertex.getPoint()); //third point of original triangle
			Point_3 p3 = new Point_3(edge.opposite.next.vertex.getPoint()); //third point of new triangle
		
			double minDistanceToP3 = lowestBoundDistance(source3D.getPoint(),p3);
			
			Matrix translationOfPointP0= translation(-p0.x,-p0.y,-p0.z);
			Point_3[] pointsArray = new Point_3[] {p0,p1,p2,p3};
			
			for(Point_3 p: pointsArray) {
				if(p.z.isNaN() || p.y.isNaN() || p.x.isNaN()) {
					throw new IllegalStateException("There is a null value in one of the points");
				}
			}
			
			for(Point_3 p: pointsArray) transform(p, translationOfPointP0);
			
			if(p0.z.isNaN()) throw new IllegalStateException("There is a null value in one of the points");
			
			double thetaRot1 = calculateAngle(p1.x,p1.y,1,0); //rotate around z until y becomes 0
			if(p1.y>0)thetaRot1*=-1; 
			Matrix rotation1ForP1 = rotation(thetaRot1, 0);
			for(Point_3 p: pointsArray) transform(p, rotation1ForP1);
			
			if(p0.z.isNaN()) throw new IllegalStateException("There is a null value in one of the points");
			
			double thetaRot2 = calculateAngle(p1.x,p1.z,1,0); // rotate around y until Z becumes 0
			if(p1.z<0)thetaRot2*=-1;// We do this because the angle should be negative if we go beyond PI
			Matrix rotation2ForP1 = rotation(thetaRot2, 2);
			for(Point_3 p: pointsArray) transform(p, rotation2ForP1);
			
			if(p0.z.isNaN()) {
				throw new IllegalStateException("There is a null value in one of the points");
			}
			
			double thetaRot3 = calculateAngle(p2.y,p2.z,1,0);
			if(p2.z>0)thetaRot3*=-1;
			Matrix rotation3ForP1 = rotation(thetaRot3, 1);
			for(Point_3 p: pointsArray) transform(p, rotation3ForP1);
			
			if(p0.z.isNaN()) throw new IllegalStateException("There is a null value in one of the points");
			
			//TODO check what happens when p2.z == 0 or something close because of EPSILON and posible error.
			
			if(p2.y<-EPSILON) {
				throw new IllegalStateException("There was an error in the calculation of the transformations");
			}
			if(p1.x<-EPSILON) throw new IllegalStateException("There was an error in the calculation of the transformations");
			if(!isAlmostZero(p0.z) || !isAlmostZero(p1.z) || !isAlmostZero(p2.z) || !isAlmostZero(p1.y)) {
//				System.out.println("p0z: " + p0 + ", p1: " + p1 + ", p2: " +p2);
				throw new IllegalStateException("There was an error in the calculation of the transformations");
			}else {
				//We are doing this to avoid problems. Maybe an error but almost nothing.
				p0.z=0.0;
				p1.z=0.0;
				p2.z=0.0;
				p1.y=0.0;
			}
			
			double distanceToLine = Math.sqrt(p3.z*p3.z + p3.y*p3.y);
			if(isAlmostZero(p2.y)) throw new IllegalStateException("there shouldn't be triangle with 3 points in the same line");
			p3.setY(distanceToLine*-1);
			p3.setZ(0);
			
			Point_3 source = calculateSourceInXWithP0in000();
			
			
			System.out.println("new Point p3(" + edge.opposite.next.vertex.getPoint() +": " + euclidianDistance(source, p3));
			
			
			Point_3 q0 = new Point_3(b0,0,0);
			Vector_3 vp0_p3 = (Vector_3) p0.minus(p3);
			Vector_3 vs_q0 = (Vector_3) source.minus(q0);
			Vector_3 vp0_s = (Vector_3) p0.minus(source);
			
			Point_3 q1 = new Point_3(b1,0,0);
			Vector_3 vp3_p1 = (Vector_3) p3.minus(p1);
			Vector_3 vs_q1 = (Vector_3) source.minus(q1);
			Vector_3 vp3_s = (Vector_3) p3.minus(source);

			// p =p0, r=p0p3, q=vs, s=vsq0, q-p = vs-p0
			Point_3[] points = new Point_3[] {p0,p1,p3,source,q0,q1};
			
			double[] t0 = new double[2];
			calculateT(t0, vp0_s, vp3_s, vs_q0, vp0_p3, vp3_p1,points, 0);
			double[] t1 = new double[2];
			calculateT(t1, vp0_s, vp3_s, vs_q1, vp0_p3, vp3_p1,points, 1);
		
			addNormalWindow(windows, t0, t1, 0, edge.opposite.next, source, p0, p3);
			addNormalWindow(windows, t0, t1, 1, edge.opposite.prev, source, p3, p1);
			
			if(isAlmostZero(b0)) {
				if(!isAlmostZero(t0[0])) {
					int index=0;
					addSaddleWindow0(windows,edge.opposite.next, index, t0, source, p0, p3, p0);
					if(isAlmostZero(t0[0]-1)) {
						index=1;
						if(!isAlmostZero(t0[index])) {
							addSaddleWindow0(windows,edge.opposite.prev, index, t0, source, p3, p1, p0);
						}
					}
				}
			}
			
			if(isAlmostZero(b1-size(edge))) {
				if(!isAlmostZero(t1[1]-1)) {
					int index=1;
					addSaddleWindow1(windows,edge.opposite.prev, index, t1, source, p3, p1, p1);
					if(isAlmostZero(t1[1])) {
						index=0;
						if(!isAlmostZero(t1[index]-1)) {
							addSaddleWindow1(windows,edge.opposite.next, index, t1, source, p0, p3, p1);
						}
					}
				}
			}
		// TODO usar TAU para decir de cual de los 2 lados venis, no a donde vas.
			return windows;
			
		}
		public Point_3 calculateSourceWithPointOnEdgeX() {
			double d = (b1-b0);
			if(isAlmostZero(d-d0-d1)) {
				//the source is on the axisX (probably new source in boundary vertex
//				System.out.println("Point aligned inside of Q0-Q1");
				return new Point_3(d0,0,0);
			}
			if(isAlmostZero(d0-d-d1)) {
				System.out.println("Point aligned outside of Q0-Q1");
				return new Point_3(d0,0,0);
			}
			if(isAlmostZero(d1-d-d0)) {
				System.out.println("Point aligned outside of Q0-Q1");
				return new Point_3(-d0,0,0);
			}
			if(isAlmostZero(d)) throw new IllegalStateException("can't calculate source if b0==b1");
			double a = (d0*d0-d1*d1+d*d)/(2*d);
			double h = Math.sqrt(d0*d0-a*a);
			if(Double.isNaN(h)) {
				throw new IllegalStateException("h is NaN");
			}
			Point_3 source = new Point_3(a,h,0);
			return source;
			//TODO use "calculatePointWithDistances
		}
		
		public Point_3 calculateSourceInXWithP0in000() {
			Point_3 source = calculateSourceWithPointOnEdgeX();
			source.x+=b0;
			return source;
		}
		
//		public Point_3 calculateSourceWithTwoPoints() {
//			Point_3[] points = new Point_3[] {edge.opposite.vertex.getPoint(), edge.vertex.getPoint()};
//			Point_3 q0 = Point_3.linearCombination(points, new Double[] {(size(edge)-b0)/size(edge),b0/size(edge)});
//			Point_3 q1 = Point_3.linearCombination(points, new Double[] {(size(edge)-b1)/size(edge),b1/size(edge)});
//			Point_3 pSameSide = edge.next.vertex.getPoint();
//			double d = (b1-b0);
//			if(isAlmostZero(d-d0-d1)) {
//				//the source is on the axisX (probably new source in boundary vertex
//				System.out.println("Source is on the edge. Do something");
//				return Point_3.linearCombination(new Point_3[] {q1, q0},new Double[] {d0,d1});
//			}
//			if(isAlmostZero(b1-b0)) throw new IllegalStateException("can't calculate source if b0==b1");
//			double a = (d0*d0-d1*d1+d*d)/(2*d);
//			double h = Math.sqrt(d0*d0-a*a);
//			if(Double.isNaN(h)) {
//				throw new IllegalStateException("h is NaN");
//			}
//			double b = d-a;
//			Point_3 p2 = Point_3.linearCombination(new Point_3[] {q1, q0},new Double[] {a,b});
//			double newX = p2.x + h* (q1.y - q0.y)/d;
//			double newY = p2.y - h* (q1.x - q0.x)/d;
//			Point_3 source1 = new Point_3(newX,newY,0);
//			double newX2 = p2.x - h* (q1.y - q0.y)/d;
//			double newY2 = p2.y + h* (q1.x - q0.x)/d;
//			Point_3 source2 = new Point_3(newX2,newY2,0);
//			
//			double valueSide = (pSameSide.x-q0.x)*(q1.y-q0.y)-(pSameSide.y-q0.y)*(q1.y-q0.y);
//			double valueSource1 = (source1.x-q0.x)*(q1.y-q0.y)-(source1.y-q0.y)*(q1.y-q0.y);
//			if(Math.signum(valueSide) == Math.signum(valueSource1)) {
//				return source1;
//			}else {
//				return source2;
//			}
//		}
		
		
		// TODO pensar en los ciclos y que onda.
		private void addSaddleWindow1(ArrayList<Window> windows, Halfedge<Point_3> currentEdge, int index, double[] t1,
				Point_3 source, Point_3 edgeOrigin, Point_3 edgeDestination, Point_3 newSource) {
			double newb0 = size(currentEdge)*t1[index];
			double newb1 = size(currentEdge);
			Point_3 newPL = Point_3.linearCombination(new Point_3[] {edgeDestination, edgeOrigin},new Double[] {t1[index],1-t1[index]});
			Point_3 newPR = edgeDestination;
			double newd0 =  euclidianDistance(newPL, newSource);
			double newd1 = euclidianDistance(newPR, newSource);
			double newSigma = sigma+euclidianDistance(newSource, source);
			Window currentWindow = new Window(newb0, newb1, newd0, newd1, newSigma, false, currentEdge, newSigma + Math.min(newd0, newd1));
			windows.add(currentWindow);
		}
		private void addSaddleWindow0(ArrayList<Window> windows, Halfedge<Point_3> currentEdge, int index, double[] t0,
				Point_3 source, Point_3 edgeOrigin, Point_3 edgeDestination, Point_3 newSource) {
			double newb0 = 0;
			double newb1 = size(currentEdge)*t0[index];
			Point_3 newPL = edgeOrigin;
			Point_3 newPR = Point_3.linearCombination(new Point_3[] {edgeDestination, edgeOrigin},new Double[] {t0[index],1-t0[index]});
			double newd0 =  euclidianDistance(newPL, newSource);
			double newd1 = euclidianDistance(newPR, newSource);
			double newSigma = sigma+euclidianDistance(newSource, source);
			Window currentWindow = new Window(newb0, newb1, newd0, newd1, newSigma, false, currentEdge, newSigma + Math.min(newd0, newd1));
			windows.add(currentWindow);
		}
		private void addNormalWindow(ArrayList<Window> windows, double[] t0, double[] t1, int index, Halfedge<Point_3> currentEdge,
				Point_3 source, Point_3 edgeOrigin, Point_3 edgeDestination) {
			if(!isAlmostZero(t1[index]-t0[index])) {
				double newb0 = size(currentEdge)*t0[index];
				double newb1 = size(currentEdge)*t1[index];
				Point_3 newPL = Point_3.linearCombination(new Point_3[] {edgeDestination, edgeOrigin},new Double[] {t0[index],1-t0[index]});
				Point_3 newPR = Point_3.linearCombination(new Point_3[] {edgeDestination, edgeOrigin},new Double[] {t1[index],1-t1[index]});
				double newd0 = euclidianDistance(newPL, source);
				double newd1 = euclidianDistance(newPR, source);
				if(Double.isNaN(newd0) || Double.isNaN(newd1)) throw new IllegalStateException("Can't have distances Nan");
				Window currentWindow = new Window(newb0, newb1, newd0, newd1, sigma, false, currentEdge, Math.min(newd0, newd1)+sigma);
				windows.add(currentWindow);
			}
			
		}
		
		public void copy(Window fw) {
			this.b0 = fw.b0;
			this.b1 = fw.b1;
			this.d0 = fw.d0;
			this.d1 = fw.d1;
			this.sigma = fw.sigma;
			this.edge = fw.edge;
			this.pqCriteria = fw.pqCriteria;
			this.tau = fw.tau;
		}
		
	}
	
	private double calculateShortestDistance(Vertex<Point_3> source, Vertex<Point_3> destination, double upperBound, Map<Halfedge<Point_3>, ArrayList<Window>> mapWindowsEdge ) {
		PriorityQueue<Window> pq = new PriorityQueue<>();
		Halfedge<Point_3> inicialHalfEdge = source.getHalfedge(); // edge that point to the vertex
		Halfedge<Point_3> currentHalfEdge = inicialHalfEdge;	
		do {
			double d0 =  size(currentHalfEdge.next); //.next is the edge that start with the vertex .prev the one that do not have the vertex
			double d1 = size(currentHalfEdge);
			Window currentW = new Window(0, size(currentHalfEdge.prev), d0, d1, 0, false, currentHalfEdge.prev, Math.min(d0, d1));
			addToQueue(pq,currentW, mapWindowsEdge);
			currentHalfEdge = currentHalfEdge.next.opposite;
		}while(inicialHalfEdge!=currentHalfEdge);
		
		int numIterationPQ = 0;
		while(!pq.isEmpty()) {
			Window currentWindow = pq.poll();
			if(currentWindow == null) throw new IllegalStateException(); // TODO delete 
			System.out.println("iteration: " + numIterationPQ++ + " " + currentWindow);
			if(isAlmostZero(currentWindow.b0)) {
				if(currentWindow.edge.opposite.vertex.equals(destination)) {
					System.out.println("destination!");
					return (currentWindow.d0+currentWindow.sigma);
				}
			}
			if(isAlmostZero(currentWindow.b1-size(currentWindow.edge))) {
				if(currentWindow.edge.vertex.equals(destination)) {
					System.out.println("destination!");
					return (currentWindow.d1+currentWindow.sigma);
				}
			}
			ArrayList<Window> windowsToAdd = currentWindow.propagate(source);
			ArrayList<Window> trueWindowsToAdd = new ArrayList();
			for(Window wToAdd: windowsToAdd) {
				Point_3 sourceOfWindowToAdd = wToAdd.calculateSourceInXWithP0in000();
				double distance=0;
//				if(sourceOfWindowToAdd.x>wToAdd.b0 && sourceOfWindowToAdd.x<wToAdd.b1) {
//					distance = Math.abs(sourceOfWindowToAdd.y);
//				}else {
//					distance = Math.min(wToAdd.d0, wToAdd.d1);
//				}
				// TODO get a more strict condition
				Point_3[] points = new Point_3[] {wToAdd.edge.opposite.vertex.getPoint(), wToAdd.edge.vertex.getPoint()};
				Point_3 q0 = Point_3.linearCombination(points, new Double[] {(size(wToAdd.edge)-wToAdd.b0)/size(wToAdd.edge),(wToAdd.b0)/size(wToAdd.edge)});
				Point_3 q1 = Point_3.linearCombination(points, new Double[] {(size(wToAdd.edge)-wToAdd.b1)/size(wToAdd.edge),wToAdd.b1/size(wToAdd.edge)});
				if(wToAdd.d0<wToAdd.d1) {
					distance = wToAdd.d0 - 2*(wToAdd.b1-wToAdd.b0)+ lowestBoundDistance(q0, destination.getPoint());
				}else {
					distance = wToAdd.d1 - 2*(wToAdd.b1-wToAdd.b0) + lowestBoundDistance(q1, destination.getPoint());
				}
				if(distance<=upperBound+EPSILON) {
					trueWindowsToAdd.add(wToAdd);
				}
			}
			propagate(trueWindowsToAdd, pq, mapWindowsEdge);
		}
		return -1;
	}

	//Intersects all windows and add the results to the pq, updating it and the map.
	private void propagate(ArrayList<Window> windowsToAdd, PriorityQueue<Window> pq,
			Map<Halfedge<Point_3>, ArrayList<Window>> mapWindowsEdge) {
		for(int listIndex=0; listIndex<windowsToAdd.size(); listIndex++) {
			Window windowInList = windowsToAdd.get(listIndex);
			if(!mapWindowsEdge.containsKey(windowInList.edge)) {
				addToQueue(pq, windowInList, mapWindowsEdge);
			}else {
				ArrayList<Window> windowsInEdge = mapWindowsEdge.get(windowInList.edge);
				for(int edgeIndex=0; edgeIndex<windowsInEdge.size(); edgeIndex++) {
					Window windowInEdge = windowsInEdge.get(edgeIndex); 
					if(windowInList.b1>windowInEdge.b0 && windowInList.b0<windowInEdge.b1) {
						Point_3 sourceList = windowInList.calculateSourceInXWithP0in000();
						Point_3 sourceEdge = windowInEdge.calculateSourceInXWithP0in000();
						// TODO check interchanging s1 and s0
						double alpha = sourceEdge.x-sourceList.x;
						double beta = windowInEdge.sigma - windowInList.sigma;
						double gama = (sourceList.x*sourceList.x+sourceList.y*sourceList.y)-(sourceEdge.x*sourceEdge.x+sourceEdge.y*sourceEdge.y)-beta*beta;
						double A = (alpha*alpha-beta*beta);
						double B = gama*alpha+2*sourceEdge.x*beta*beta;
						double C = 1/4*gama*gama-(sourceEdge.x*sourceEdge.x+sourceEdge.y*sourceEdge.y)*beta*beta;
						
						double[] roots = new double[2];
						int numRoots = caculateRoots(roots, A, B, C);
						boolean isThereARoot = false;
						double rootInTheSegment = -1;
						
						double leftBound = windowInEdge.b0>windowInList.b0?windowInEdge.b0:windowInList.b0;
						double rightBound = windowInEdge.b1<windowInList.b1?windowInEdge.b1:windowInList.b1;
						if(numRoots>0) {
							boolean isFirstRootIn = isRootInIntersection(roots[0],leftBound,rightBound);
							boolean isSecondRootIn = false;
							if(numRoots==2) {
								isSecondRootIn = isRootInIntersection(roots[1],leftBound,rightBound);
							}
							if(isFirstRootIn && isSecondRootIn) {
								throw new IllegalStateException("There are two roots in the interval");
								//We could solve this easy just do 3 windows and two outwindows but it is another case. Not dificult... but is it normal?
							}
							if(isFirstRootIn) {
								isThereARoot=true;
								rootInTheSegment=roots[0];
							}else if(isSecondRootIn) {
								isThereARoot=true;
								rootInTheSegment=roots[1];
							}
						}
						ArrayList<Window> newWindowsListAfterIntersection = new ArrayList<>();
						ArrayList<Window> newWindowsEdgeAfterIntersection = new ArrayList<>();
						if(!isThereARoot) { //All the segment correspond to one of the windows.
							double middleOfSegment = (rightBound+leftBound)/2;
							Point_3 middlePoint = new Point_3(middleOfSegment,0,0);
							if(euclidianDistance(middlePoint, sourceList) + windowInList.sigma<euclidianDistance(middlePoint, sourceEdge) + windowInEdge.sigma) {
								//windowListIsTheSmallerDistance Strict smaller so we don't have loops
								calculateNewRightWindow(newWindowsEdgeAfterIntersection, windowInEdge, windowInList, sourceEdge);
								calculateNewLeftWindow(newWindowsEdgeAfterIntersection, windowInEdge, windowInList, sourceEdge);
								newWindowsListAfterIntersection.add(windowInList);
							}else { //windowEdgeIsTheSmallerDistance
								calculateNewRightWindow(newWindowsListAfterIntersection, windowInList, windowInEdge, sourceList);
								calculateNewLeftWindow(newWindowsListAfterIntersection, windowInList, windowInEdge, sourceList);
							}
						}else {
							Point_3 pointToLeftOfRoot = new Point_3((rootInTheSegment+leftBound)/2,0,0);
							Point_3 pointToRightOfRoot = new Point_3((rootInTheSegment+rightBound)/2,0,0);
							if(euclidianDistance(pointToLeftOfRoot, sourceList) + windowInList.sigma<euclidianDistance(pointToLeftOfRoot, sourceEdge) + windowInEdge.sigma) {
								//windowLISTIsTheSmallerDistance on the LEFT of the root
								calculateNewWindows(newWindowsListAfterIntersection, newWindowsEdgeAfterIntersection, windowInList, windowInEdge, sourceList, sourceEdge, rootInTheSegment);
							}else {
								//windowEDGEIsTheSmallerDistance on the LEFT of the root
								calculateNewWindows(newWindowsEdgeAfterIntersection, newWindowsListAfterIntersection, windowInEdge, windowInList, sourceEdge, sourceList, rootInTheSegment);
							}
						}
						//TODO add windows. check que el iterator no explote.
						if(newWindowsListAfterIntersection.size()>0) {
							Window fw = newWindowsListAfterIntersection.get(0);
							windowInList.copy(fw);
							if(windowInList.b1-windowInList.b0 <EPSILON) {
								throw new IllegalStateException("We are adding window super small");
							}
							for(int i=1; i<newWindowsListAfterIntersection.size(); i++) {
								if(newWindowsListAfterIntersection.get(i).b1-newWindowsListAfterIntersection.get(i).b0 <EPSILON) {
									throw new IllegalStateException("We are adding window super small");
								}
								windowsToAdd.add(newWindowsListAfterIntersection.get(i));
							}
						}else {
							windowInList.b1=windowInList.b0;
							break; // TODO do not use break;
						}
						//TODO do not repeat code
						boolean wasInPQ = pq.remove(windowInEdge);
						if(newWindowsEdgeAfterIntersection.size()>0) {
							Window fw = newWindowsEdgeAfterIntersection.get(0);
							windowInEdge.copy(fw);
							if(wasInPQ) {
								if(windowInEdge.b1-windowInEdge.b0 <EPSILON) {
									throw new IllegalStateException("We are adding to PQ something super small");
								}
								pq.add(windowInEdge);
								
							}
							for(int i=1; i<newWindowsEdgeAfterIntersection.size(); i++) {
								Window w = newWindowsEdgeAfterIntersection.get(i);
								if(w.b1-w.b0 <EPSILON) {
									throw new IllegalStateException("We are adding to PQ something super small");
								}
								windowsInEdge.add(w);
								if(wasInPQ) pq.add(w);
							}
						}else {
							windowsInEdge.remove(windowInEdge);
							edgeIndex--;
						}
					}
				}
				if(!isAlmostZero(windowInList.b0-windowInList.b1)) {
					pq.add(windowInList);
					windowsInEdge.add(windowInList);
				}
			}
			ArrayList<Window> windowsInArrayEdge = mapWindowsEdge.get(windowInList.edge);
			for(Window w: windowsInArrayEdge) {
				for(Window w2 : windowsInArrayEdge) {
					if(w!=w2) {
						if(w.b1>w2.b0 && w.b0<w2.b1) {
							throw new IllegalStateException("Error calculation intersections");
							// TODO delete this for optimization.
						}
					}
				}
			}
		}
		
	}

	private void calculateNewWindows(ArrayList<Window> newWindowsLeftBetter,
			ArrayList<Window> newWindowsRightBetter, Window windowLeftBetter, Window windowRightBetter,
			Point_3 sourceLeftBetter, Point_3 sourceRightBetter, double root) {
		// TODO Auto-generated method stub
		Point_3 rootPoint = new Point_3(root,0,0);
		double newD1 = euclidianDistance(rootPoint, sourceLeftBetter);
		Window leftWindow = new Window(windowLeftBetter.b0, root, windowLeftBetter.d0, newD1, windowLeftBetter.sigma, windowLeftBetter.tau, windowLeftBetter.edge, Math.min(windowLeftBetter.d0, newD1)+windowLeftBetter.sigma);
		if(leftWindow.b1-leftWindow.b0 <EPSILON) {
			throw new IllegalStateException();
		}
		newWindowsLeftBetter.add(leftWindow);
		calculateNewLeftWindow(newWindowsRightBetter, windowRightBetter, windowLeftBetter, sourceRightBetter);
		double newD0 = euclidianDistance(rootPoint, sourceRightBetter);
		
		Point_3 rightestPoint = new Point_3(windowRightBetter.b1,0,0);
		double rightestPointDistance = euclidianDistance(rightestPoint, sourceRightBetter);
		if(!isAlmostZero(rightestPointDistance-windowRightBetter.d1)) {
			throw new IllegalStateException();
		}
		
		Window rightWindow = new Window(root, windowRightBetter.b1, newD0, windowRightBetter.d1, windowRightBetter.sigma, windowRightBetter.tau, windowRightBetter.edge, Math.min(windowRightBetter.d1, newD0)+windowRightBetter.sigma);
		if(rightWindow.b1-rightWindow.b0 <EPSILON) {
			throw new IllegalStateException();
		}
		newWindowsRightBetter.add(rightWindow);
		calculateNewRightWindow(newWindowsLeftBetter, windowLeftBetter, windowRightBetter, sourceLeftBetter);
	}

	private void calculateNewLeftWindow(ArrayList<Window> newWindowsAfterIntersection, Window windowWorse,
			Window windowBetter, Point_3 sourceWorse) {
		Point_3 leftPointBetter = new Point_3(windowBetter.b0,0,0);
		Point_3 leftPointWorse = new Point_3(windowWorse.b0,0,0);
		if(!isAlmostZero(windowWorse.d0-euclidianDistance(leftPointWorse, sourceWorse))) {
			throw new IllegalStateException("error distances for extra window calculation");
		}
		if(windowWorse.b0<windowBetter.b0 && !isAlmostZero(windowWorse.b0-windowBetter.b0)) {
			Window newWindow = new Window(windowWorse.b0, windowBetter.b0, windowWorse.d0, euclidianDistance(leftPointBetter, sourceWorse), windowWorse.sigma, windowWorse.tau, windowWorse.edge, Math.min(windowWorse.d0, euclidianDistance(leftPointBetter, sourceWorse))+windowWorse.sigma);
			if(newWindow.b1-newWindow.b0 <EPSILON) {
				throw new IllegalStateException();
			}
			newWindowsAfterIntersection.add(newWindow);
		}
	}
	
	private void calculateNewRightWindow(ArrayList<Window> newWindowsAfterIntersection, Window windowWorse,
			Window windowBetter, Point_3 sourceWorse) {
		Point_3 rightPointBetter = new Point_3(windowBetter.b1,0,0);
		Point_3 rightPointWorse = new Point_3(windowWorse.b1,0,0);
		if(!isAlmostZero(windowWorse.d1-euclidianDistance(rightPointWorse, sourceWorse))) throw new IllegalStateException("error distances for extra window calculation");
		if(windowWorse.b1>windowBetter.b1 && !isAlmostZero(windowWorse.b1-windowBetter.b1)) {
			Window newWindow = new Window(windowBetter.b1, windowWorse.b1, euclidianDistance(rightPointBetter, sourceWorse), windowWorse.d1, windowWorse.sigma, windowWorse.tau, windowWorse.edge, Math.min(euclidianDistance(rightPointBetter, sourceWorse), windowWorse.d1)+windowWorse.sigma);
			if(newWindow.b1-newWindow.b0 <EPSILON) {
				throw new IllegalStateException();
			}
			newWindowsAfterIntersection.add(newWindow);
		}
	}

	private boolean isRootInIntersection(double root, double left, double right) {
		return (!isAlmostZero(left-root) && !isAlmostZero(right-root) && (root>left && root<right));
	}

	private int caculateRoots(double[] roots, double a, double b, double c) {
		  double d = b * b - 4 * a * c;
	        if(d > 0) {
	            roots[0] = ( - b + Math.sqrt(d))/(2*a);
	            roots[1] = (-b - Math.sqrt(d))/(2*a);
	            return 2;
	        }
	        else if(d == 0){
	            roots[0] = (-b+Math.sqrt(d))/(2*a);
	            return 1;
	        }
	        return 0;
	}

	private void addToQueue(PriorityQueue<Window> pq, Window currentW,
			Map<Halfedge<Point_3>, ArrayList<Window>> mapWindowsEdge) {
		
		if(!mapWindowsEdge.containsKey(currentW.edge)) {
			mapWindowsEdge.put(currentW.edge, new ArrayList<>());
		}
		mapWindowsEdge.get(currentW.edge).add(currentW);
		pq.add(currentW);
	}

	private static Matrix rotation(double theta, int offset) { //offset= 0 for z, 1 for x and 2 for y
		double[][] array = new double[4][4];
		// TODO delete this is java not C
		for(int i=0; i<4; i++) {
			for(int j=0; j<4; j++) {
				array[i][j] = 0;
			}
		}
		array[3][3] = 1;
		array[(2+offset)%3][(2+offset)%3]=1;
		array[(1+offset)%3][(1+offset)%3]=Math.cos(theta);
		array[(0+offset)%3][(0+offset)%3]=Math.cos(theta);
		array[(0+offset)%3][(1+offset)%3]=-Math.sin(theta);
		array[(1+offset)%3][(0+offset)%3]=Math.sin(theta);
		return new Matrix(array);
	}

	private static double calculateAngle(double x, double y, double i, double j) {
		if(isAlmostZero(y)) {
			if(x>=0) return 0;
			return Math.PI;
		}
		double a = Math.sqrt(x*x+y*y);
		double c = Math.sqrt((x-i)*(x-i) + (y-j)*(y-j));
		return Math.acos((1+a*a-c*c)/(2*a));
	}

	private static Matrix translation(Double x, Double y, Double z) {
		double[][] array = {{1.,0., 0.,x},{0.,1.,0., y},{0.,0.,1.,z},{0.,0.,0.,1}};
		Matrix t = new Matrix(array);
		return t;
	}
	/**
	 * Apply the transformation to point p (having homogeneous coordinates)
	 */
	public static void transform(Point_3 p, Matrix m) {
		double x=p.getX().doubleValue();
		double y=p.getY().doubleValue();
		double z=p.getZ().doubleValue();
		double[][] array = {{x}, {y}, {z}, {1}}; 
		Matrix v=new Matrix(array); // the vector
		
		Matrix result=m.times(v);
		p.setX(result.get(0, 0));
		p.setY(result.get(1, 0));
		p.setZ(result.get(2, 0));
		//return new Point_3(result.get(0, 0), result.get(1, 0), result.get(2, 0));
	}


	// TODO put this in the information of the halfedge to make it faster... maybe?
	private static double size(Halfedge<Point_3> edge) {
		return euclidianDistance(edge.vertex.getPoint(), edge.opposite.vertex.getPoint());
	}
	
	private static double lowestBoundDistance(Point_3 source, Point_3 destination) {
		return euclidianDistance(source, destination);
		// TODO step 2
	}
	
	// DIJKSTRA
	private static class DijkstraNode implements Comparable<DijkstraNode>{
		Vertex<Point_3> v;
		double distance;
		int depth; // TODO remove
		public DijkstraNode(Vertex<Point_3> v, double distance, int depth) {
			super();
			this.v = v;
			this.distance = distance;
			this.depth = depth;
		}
		@Override
		public int compareTo(DijkstraNode dn) {
			if(dn == null) throw new IllegalArgumentException();
			double diff = distance-dn.distance;
			if(diff>0) return 1;
			if(diff<0) return -1;
			return 0;
		}
		@Override
		public int hashCode() {
			final int prime = 31;
			int result = (int)distance;
			result = prime * result + ((v == null) ? 0 : v.hashCode());
			return result;
		}
		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			DijkstraNode other = (DijkstraNode) obj;
			if (v == null) {
				if (other.v != null)
					return false;
			} else if (v!=other.v)
				return false;
			return true;
		}
	}

	public double bidirectionalDijkstra(Vertex<Point_3> source, Vertex<Point_3> destination) {
		//Source Dijkstra variables
		PriorityQueue<DijkstraNode> pqSource = new PriorityQueue<>();
		pqSource.add(new DijkstraNode(source, 0, 0));
		Map<Vertex<Point_3>, Double> mapOfNodesSource = new HashMap<Vertex<Point_3>, Double>();
		
		//Destination Dijkstra variables
		PriorityQueue<DijkstraNode> pqDestination = new PriorityQueue<>();
		pqDestination.add(new DijkstraNode(destination, 0, 0));
		Map<Vertex<Point_3>, Double> mapOfNodesDestination = new HashMap<Vertex<Point_3>, Double>();
		
		boolean isSourceTurn = true;
		double posibleMinDistance = -1;
		while(true) {
			if(isSourceTurn) {
				posibleMinDistance = doDijkstraStep(pqSource, mapOfNodesSource, mapOfNodesDestination);	
			}else {
				posibleMinDistance = doDijkstraStep(pqDestination, mapOfNodesDestination, mapOfNodesSource);	
			}
			if(posibleMinDistance>=0) {
				return posibleMinDistance;
			}
			isSourceTurn = !isSourceTurn;
		}
	}

	private double doDijkstraStep(PriorityQueue<DijkstraNode> pqSource, Map<Vertex<Point_3>, Double> mapOfNodesSource,
			Map<Vertex<Point_3>, Double> mapOfNodesDestination) {
			if(!pqSource.isEmpty()) {
				DijkstraNode currentNode = pqSource.poll();
				if(mapOfNodesDestination.containsKey(currentNode.v)) {
					System.out.println("depth of 1 of the dijkstra is: " + currentNode.depth);
					return currentNode.distance + mapOfNodesDestination.get(currentNode.v);
				}
				if(!mapOfNodesSource.containsKey(currentNode.v)){
					mapOfNodesSource.put(currentNode.v, currentNode.distance);
					addNeighbors(currentNode, pqSource, mapOfNodesSource);
				}
			}else {
				throw new IllegalStateException("The two vertices are not connected");
			}
			return -1;
	}

	private void addNeighbors(DijkstraNode currentNode, PriorityQueue<DijkstraNode> pq, Map<Vertex<Point_3>, Double> mapOfNodesSource) {
		Vertex<Point_3> inicialV = currentNode.v;
		double inicialDistance = currentNode.distance;
		Halfedge<Point_3> inicialHalfEdge = inicialV.getHalfedge();
		Halfedge<Point_3> currentHalfEdge = inicialHalfEdge;	
		do {
			Vertex<Point_3> currentV = currentHalfEdge.opposite.vertex;
			if(!mapOfNodesSource.containsKey(currentV)) {
				double currentDistance = inicialDistance+ euclidianDistance(inicialV.getPoint(), currentV.getPoint());
				pq.add(new DijkstraNode(currentV, currentDistance, currentNode.depth+1));
			}
			currentHalfEdge = currentHalfEdge.next.opposite;
		}while(inicialHalfEdge!=currentHalfEdge);
	}
	//END OF DIJKSTRA

	private static double euclidianDistance(Point_3 inicialP, Point_3 currentP) {
		double value =  Math.sqrt(inicialP.squareDistance(currentP).doubleValue());
		if(Double.isNaN(value)) {
			throw new IllegalArgumentException("Can't have distance nan");
		}
		return value;
	}

	public static boolean isAlmostZero(Double d) {
		return d<EPSILON && d>-EPSILON;
	}
//	calculateT(t, vp1_p, vp2_p, vp_s, vp1_p2, vp2_p0,points, 0);
	
//	public static void calculateTBackPropagation(double[] t, Vector_3 vDiffOrigin, Vector_3 vDiffOrigin2, Vector_3 vLineSource, Vector_3 vLineLeft,
//			Vector_3 vLineRight, Point_3[] points, int index) { //Points:p0-p1-p3-source-q0-q1
//		double auxCross = cross(vLineLeft, vLineSource);
//		if(isAlmostZero(auxCross)) {
//			System.out.println("Paralel lines");
//			if(isAlmostZero(cross(vDiffOrigin,vLineSource))) {
//				t[0] = 1;
//			}else {
//				t[0] = 2;
//			}
//		}else {
//			t[0] = cross(vDiffOrigin,vLineSource)/auxCross;
//		}
//		double auxCross2 = cross(vLineRight, vLineSource);
//		if(isAlmostZero(auxCross2)) {
//			t[1] = 2;
//		}else {
//			t[1] = cross(vDiffOrigin2,vLineSource)/auxCross2;
//		}
//		//if it cuts both lines:
//		if((t[0]>=-EPSILON && t[0]<=1+EPSILON) && (t[1]>=-EPSILON && t[1]<=1+EPSILON)) {
//			if(!isAlmostZero(t[1])) {
//				if(isAlmostZero(t[0])) t[0] = 1;
//				else {
//					if(isAlmostZero(t[1]-1)) t[1] = 0;
//				}
//			}
//		}else {
//			if(t[0]<-EPSILON || t[0]>1+EPSILON) {
//				t[0]=1;
//			}
//			if(t[1]<-EPSILON || t[1]>1+EPSILON) {
//				t[1] = 0;
//			}
//		}
//		if((t[0]<-EPSILON || t[0]>1+EPSILON) && (t[1]<-EPSILON || t[1]>1+EPSILON)) {
//			throw new IllegalArgumentException("needs to cut one of the two lines");
//		}
//		if((t[0]>=EPSILON && t[0]<=1-EPSILON) && (t[1]>=EPSILON && t[1]<=1-EPSILON)) {
//			throw new IllegalStateException("can't cut both lines in the middle");
//		}
//	}
	public static void calculateT(double[] t, Vector_3 vDiffOrigin, Vector_3 vDiffOrigin2, Vector_3 vLineSource, Vector_3 vLineLeft,
			Vector_3 vLineRight, Point_3[] points, int index) { //Points:p0-p1-p3-source-q0-q1
		// Fix for when q0==p0 || q1==p1
//		if(points.length!=6) throw new IllegalStateException("Point list needs to be of lenght 6");
//		if(isAlmostZero(points[0+index].x-points[4+index].x)) { //MAL no es siempre verdad y con esto solo podes mirar ventanas enteras:(
//			t[0]=index;
//			t[1]=index;
//			return;
//		}
		boolean isParallel = false;
		double auxCross = cross(vLineLeft, vLineSource);
		if(isAlmostZero(auxCross)) {
			System.out.println("Paralel lines 1");
			t[0] = 2;
			isParallel = true;
		}else {
			t[0] = cross(vDiffOrigin,vLineSource)/auxCross;
		}
		double auxCross2 = cross(vLineRight, vLineSource);
		if(isAlmostZero(auxCross2)) {
			System.out.println("Paralel lines 2");
			t[1] = -1;
			if(isParallel) {
				throw new IllegalStateException("Can't be paralel to both lines");
			}
		}else {
			t[1] = cross(vDiffOrigin2,vLineSource)/auxCross2;
		}
		//if it cuts both lines:
		if((t[0]>=-EPSILON && t[0]<=1+EPSILON) && (t[1]>=-EPSILON && t[1]<=1+EPSILON)) {
			if(isAlmostZero(t[0]) && isAlmostZero(t[1]-1)) {
				//Points:p0-p1-p3-source-q0-q1
//				throw new IllegalStateException("border case. Add something in case it happens");
				t[0] = 1;
				t[1] = 0;
				System.out.println("border case. Can be wrong");
			}else {
				if(isAlmostZero(t[0])) {
					t[0]=1;
				}
				if(isAlmostZero(t[1]-1)) {
					t[1]=0;
				}
			}
		}else {
			if((t[0]<-EPSILON || t[0]>1+EPSILON) && (t[1]<-EPSILON || t[1]>1+EPSILON)) {
				throw new IllegalArgumentException("needs to cut one of the two lines");
			}
			if(t[0]<-EPSILON || t[0]>1+EPSILON) {
				//No corta la primera linea.
				t[0] = 1;
			}
			if(t[1]<-EPSILON || t[1]>1+EPSILON) {
				//No corta la segunda
				t[1] = 0;
			}
		}
		if((t[0]<-EPSILON || t[0]>1+EPSILON) && (t[1]<-EPSILON || t[1]>1+EPSILON)) {
			throw new IllegalArgumentException("needs to cut one of the two lines");
		}
		if((t[0]>=EPSILON && t[0]<=1-EPSILON) && (t[1]>=EPSILON && t[1]<=1-EPSILON)) {
			throw new IllegalStateException("can't cut both lines in the middle");
		}
	}
	public static double cross(Vector_3 v1, Vector_3 v2) {
		return v1.x*v2.y - v1.y*v2.x;
	}
	

}
