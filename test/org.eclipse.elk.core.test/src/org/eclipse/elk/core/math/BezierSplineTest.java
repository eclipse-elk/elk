package org.eclipse.elk.core.math;

import static org.junit.Assert.assertEquals;

import org.junit.Test;

public class BezierSplineTest {

    private BezierSpline bs;
    private KVector startPnt;
    private KVector fstCtrPnt;
    private KVector sndCtrPnt;
    private KVector endPnt;

    /**
     * initialize a new BezierSpline with Orizantal ranged KVectors
     */
    private void initBezierCurvehOrizontal() {
        startPnt = new KVector(0, 50);
        fstCtrPnt = new KVector(10, 50);
        sndCtrPnt = new KVector(20, 50);
        endPnt = new KVector(30, 50);
        bs = new BezierSpline();
        bs.addCurve(startPnt, fstCtrPnt, sndCtrPnt, endPnt);
    }

    /**
     * initialize a new BezierSpline with Vertical ranged KVectors
     */
    private void initBezierCurveVertical() {
        startPnt = new KVector(50, 0);
        fstCtrPnt = new KVector(50, 20);
        sndCtrPnt = new KVector(50, 40);
        endPnt = new KVector(50, 50);
        bs = new BezierSpline();
        bs.addCurve(startPnt, fstCtrPnt, sndCtrPnt, endPnt);
    }

    /**
     * 
     * Test getStartPoint from BezierSpline class.
     * 
     */
    @Test
    public void testGetStartPoint() {
        this.initBezierCurveVertical();
        KVector v = bs.getStartPoint();
        assertEquals(startPnt.x, v.x, 0);
        assertEquals(startPnt.y, v.y, 0);
        this.initBezierCurvehOrizontal();
        v = bs.getStartPoint();
        assertEquals(startPnt.x, v.x, 0);
        assertEquals(startPnt.y, v.y, 0);
    }

    /**
     * 
     * Test getEndPoint from BezierSpline class.
     * 
     */
    @Test
    public void testGetEndPoint() {
        this.initBezierCurveVertical();
        KVector v = bs.getEndPoint();
        assertEquals(endPnt.x, v.x, 0);
        assertEquals(endPnt.y, v.y, 0);
        this.initBezierCurvehOrizontal();
        v = bs.getEndPoint();
        assertEquals(endPnt.x, v.x, 0);
        assertEquals(endPnt.y, v.y, 0);

    }

    /**
     * 
     * Test getInnerPoints from BezierSpline class.
     * 
     */
    @Test
    public void testGetInnerPoints() {
        // test if all InnerPoints have the same x when a Vertical Curve
        this.initBezierCurveVertical();
        KVector[] vectors = bs.getInnerPoints();
        for (KVector v : vectors) {
            assertEquals(startPnt.x, v.x, 0.00000001);
        }
        
        // test if all InnerPoints have the same y when a Orizontal Curve
        this.initBezierCurvehOrizontal();
        vectors = bs.getInnerPoints();
        for (KVector v : vectors) {
            assertEquals(startPnt.y, v.y, 0.00000001);
        }
    }

    /**
     * 
     * Test getBasePoints from BezierSpline class.
     * 
     */
    @Test
    public void testGetBasePoints() {
        // test if all BasePoints have the same x when a Vertical Curve
        this.initBezierCurveVertical();
        KVector[] vectors = bs.getBasePoints();
        for (KVector v : vectors) {
            assertEquals(startPnt.x, v.x, 0.00000001);
        }
        
        // test if all InnerPoints have the same y when a Orizontal Curve
        this.initBezierCurvehOrizontal();
        vectors = bs.getBasePoints();
        for (KVector v : vectors) {
            assertEquals(startPnt.y, v.y, 0.00000001);
        }
    }

    /**
     * 
     * Test getPolylineApprx from BezierSpline class.
     * 
     */
    @Test
    public void testGetPolylineApprx() {
        // test if all PolylineApprx have the same x when a Vertical Curve
        this.initBezierCurveVertical();
        KVector[] vectors = bs.getPolylineApprx(50);
        for (KVector v : vectors) {
            assertEquals(startPnt.x, v.x, 0.00000001);
        }
        
        // test if all PolylineApprx have the same y when a Orizontal Curve
        this.initBezierCurvehOrizontal();
        vectors = bs.getPolylineApprx(50);
        for (KVector v : vectors) {
            assertEquals(startPnt.y, v.y, 0.00000001);
        }
    }

}
