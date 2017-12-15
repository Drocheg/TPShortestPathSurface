import processing.core.*;

import java.util.ArrayList;
import java.util.concurrent.ThreadLocalRandom;

import Jcg.geometry.*;
import Jcg.polyhedron.*;

/**
 * A simple 3d viewer for visualizing surface meshes (based on Processing)
 * 
 * @author Luca Castelli Aleardi (INF555, 2012)
 *
 */
public class MeshViewer extends PApplet {

    SurfaceMesh mesh; // 3d surface mesh
    int renderType=0; // choice of type of rendering
    int renderModes=3; // number of rendering modes
    int renderPathType=0;
    int renderPathModes=3;
    ShortestPathCalculator spc;
    Vertex<Point_3> s;
    Vertex<Point_3> d;
    ArrayList<Point_3> shortesPath;

    //String filename="OFF/high_genus.off";
//    String filename="OFF/sphere.off";
    //	String filename="OFF/cube.off";
    //	String filename="OFF/torus_33.off";
    	String filename="OFF/tore.off";
//    	String filename="OFF/cow.off";
    //String filename="OFF/letter_a.off";
    //	String filename="OFF/star.off";
    //String filename="OFF/tri_triceratops.off";
    ArrayList<String> filenames;
    int fileIndex = 3;

    public void setup() {
        filenames = new ArrayList();
        filenames.add("OFF/sphere.off");
        filenames.add("OFF/cube.off");
        filenames.add("OFF/star.off");
        filenames.add("OFF/tore.off");
        fileIndex = 0;
        size(800,600,P3D);
        ArcBall arcball = new ArcBall(this); // TODO delete?

        this.mesh=new SurfaceMesh(this, filename);
        spc = new ExactGeodesics(this.mesh.polyhedron3D);
        //ms.polyhedron3D.isValid(false);
    }

    public void draw() {
        background(0);
        //this.lights();
        directionalLight(101, 204, 255, -1, 0, 0);
        directionalLight(51, 102, 126, 0, -1, 0);
        directionalLight(51, 102, 126, 0, 0, -1);
        directionalLight(102, 50, 126, 1, 0, 0);
        directionalLight(51, 50, 102, 0, 1, 0);
        directionalLight(51, 50, 102, 0, 0, 1);

        translate(width/2.f,height/2.f,-1*height/2.f);
        this.strokeWeight(1);
        stroke(150,150,150);

        this.mesh.draw(renderType, renderPathType, s, d, shortesPath);
    }

    public void keyPressed(){
        switch(key) {
        case('d'):case('D'): this.randomDistance(); break;
        case('r'):this.renderType=(this.renderType+1)%this.renderModes; break;
        case('p'):this.renderPathType=(this.renderPathType+1)%this.renderPathModes; break;
        case('f'):{
                filename = filenames.get((fileIndex));
                fileIndex = (fileIndex+1)%filenames.size();
                this.mesh=new SurfaceMesh(this, filename);
                spc = new ExactGeodesics(this.mesh.polyhedron3D);
            }
        }
    }

    private void randomDistance() {
        ArrayList<Vertex<Point_3>> vertices = this.mesh.polyhedron3D.vertices;
        int randomNum = ThreadLocalRandom.current().nextInt(0, vertices.size());
        int randomNum2;
        do {
            randomNum2 = ThreadLocalRandom.current().nextInt(0, vertices.size());
        }while(randomNum2==randomNum);


        //TODO delete:

        // BUG SPHERE :  154 - 59
        // SPHERE: 87 - 113
        // 5 - 74
        
        // 130 - 82
        // 33 - 31

        //			randomNum = 16;
        //			randomNum2 = 40; 

        //			randomNum = 87;
        //			randomNum2 = 113;
//        randomNum = 5;
//        randomNum2 = 74;
        
        //157and 2 when not stopping Can't be paralel to both lines
        //127and 62 there was no edge that touched the destination
        //13and 97 It is not respecting the bounds
//        randomNum = 102;
//        randomNum2 = 55;
//        using 138and 26
//        It is not respecting the bounds
        //149and 80 hace una linea recta hasta un vertice y clava un angula 90° y llega al resultado
        //using 118and 90 facha la zeta al final
        
        //TORE:
//        //28and 15
//        randomNum = 28;
//        randomNum2 = 15;
//        
//        randomNum = 41;
//        randomNum2 = 5;
        //40and 26 tore  Can't find two windows in same edge with same point
        Vertex<Point_3> vSource = vertices.get(randomNum);
        Vertex<Point_3> vDestination =  vertices.get(randomNum2);
        

        System.out.println("Random numbers being used: " + randomNum + " - " + randomNum2);
        s = vSource;
        d = vDestination;
        renderType = 3;
        renderPathType = 2;
        ArrayList<Point_3> listOfPoints = new ArrayList<>();
        listOfPoints.add(vSource.getPoint());
        listOfPoints.add(vDestination.getPoint());
        try {
            listOfPoints = spc.calculatesShortestPath(vSource, vDestination);
        } catch (Exception e) {
            System.out.println("using " + randomNum + "and " + randomNum2);
            System.err.println(e.getMessage());
            System.err.println(e.getStackTrace());
            throw e;
        }
        shortesPath = listOfPoints;
        System.out.println("using " + randomNum + "and " + randomNum2);
    }

    // TODO delete
    public void subdivide() {

        this.mesh.updateScaleFactor();
        this.mesh.polyhedron3D.isValid(false);
    }

    /**
     * For running the PApplet as Java application
     */
    public static void main(String args[]) {
        //PApplet pa=new MeshViewer();
        //pa.setSize(400, 400);
        PApplet.main(new String[] { "MeshViewer" });
    }

}