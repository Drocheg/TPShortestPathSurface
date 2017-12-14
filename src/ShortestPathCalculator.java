import java.util.ArrayList;

import Jcg.geometry.Point_3;
import Jcg.polyhedron.Polyhedron_3;
import Jcg.polyhedron.Vertex;

/**
 * Calculates the shortest path between vertices.
 * Only in a triangle mesh with or without boundaries.
 * @author Droche
 *
 */
public abstract class ShortestPathCalculator {
	public Polyhedron_3<Point_3> polyhedron3D;
	
	public ShortestPathCalculator(Polyhedron_3<Point_3> polyhedron3D) {
		this.polyhedron3D=polyhedron3D;
	}
	/**
	 * The main method calculating the shortest path between two vertices.
	 */
	public abstract ArrayList<Point_3> calculatesShortestPath(Vertex<Point_3> source, Vertex<Point_3> destination);

	/**
	 * The method that calculates the shortestPath but only returns shortestDistance
	 */
	// public abstract double calculatesShortestDistance(Vertex<Point_3> source, Vertex<Point_3> destination);
	// TODO delete
}
