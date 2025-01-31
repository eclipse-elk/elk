package org.eclipse.elk.core.math;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class KVectorTest {

    /**
     * Test equals from KVector class.
     */
    @Test
    public void testEquals() {
        // init 2 KVectors
        KVector kvector1;
        KVector kvector2;
        kvector1 = new KVector(10, 10);
        kvector2 = new KVector();
        kvector2.x = 10;
        kvector2.y = 10;
        
        // test if kvector1 equals to kvector2
        assertTrue(kvector1.equals(kvector2));
        
        // init obj
        kvector1 = new KVector();
        Object obj = new Object();
        
        // test if kvector1 not equals obj
        assertFalse(kvector1.equals(obj));
    }

    /**
     * Test add and sub from KVector class.
     */
    @Test
    public void testAddAndSub() {
        // init 2 KVectors
        KVector kvector1;
        KVector kvector2;
        kvector1 = new KVector(12, 70);
        kvector2 = new KVector(15, 17);
        
        // test Add and Sub functions
        // adding and subtracting the kvector2 to kvector1 = kvector1
        assertTrue(kvector1.equals(kvector1.add(kvector2).sub(kvector2)));
    }

    /**
     * Test scale from KVector class.
     */
    @Test
    public void testScale() {
        // init KVectors
        KVector a = new KVector(12, 70);
        KVector b = new KVector(12, 70);

        KVector a_temp = new KVector(12, 70);

        a = a.add(a_temp).add(a_temp);
        b = b.scale(3);
        
        // test Scale functions
        assertTrue(a.equals(b));
    }

    /**
     * Test translate from KVector class.
     */
    @Test
    public void testTranslate() {
        // init KVector
        KVector v = new KVector(10, 30);
        KVector b = new KVector(50, 50);
        assertTrue(b.equals(v.add(40, 20)));
    }

    /**
     * Test normalize from KVector class.
     */
    @Test
    public void testNormalize() {
        // init a KVectors
        KVector v = new KVector(2, 0);
        KVector n = new KVector(1, 0);
        assertTrue(n.equals(v.normalize()));
        v = new KVector(0, 2);
        n = new KVector(0, 1);
        assertTrue(n.equals(v.normalize()));
    }

    /**
     * Test toDegrees from KVector class.
     */
    @Test
    public void testToDegrees() {
        KVector v;
        
        // test 0 degrees
        v = new KVector(10, 0);
        assertEquals(0, v.toDegrees(), 0.00001);
        
        // test 45 degrees
        v = new KVector(10, 10);
        assertEquals(45, v.toDegrees(), 0.00001);
        
        // test 90 degrees
        v = new KVector(0, 10);
        assertEquals(90, v.toDegrees(), 0.00001);
        
        // test 135 degrees
        v = new KVector(-10, 10);
        assertEquals(135, v.toDegrees(), 0.00001);
        
        // test 180 degrees
        v = new KVector(-10, 0);
        assertEquals(180, v.toDegrees(), 0.00001);
        
        // test 225 degrees
        v = new KVector(-10, -10);
        assertEquals(225, v.toDegrees(), 0.00001);
        
        // test 270 degrees
        v = new KVector(0, -10);
        assertEquals(270, v.toDegrees(), 0.00001);
        
        // test 315 degrees
        v = new KVector(10, -10);
        assertEquals(315, v.toDegrees(), 0.00001);
    }

    /**
     * Test distance from KVector class.
     */
    @Test
    public void testDistance() {
        // init 2 KVectors
        KVector v1 = new KVector(5, 50);
        KVector v2 = new KVector(5, 50);
        assertEquals(0, v1.distance(v2), 0);
        v1 = new KVector(0, 20);
        v2 = new KVector(0, 50);
        assertEquals(30, v1.distance(v2), 0);
    }

    /**
     * Test parse from KVector class.
     */
    @Test
    public void testParse() {
        // init KVector
        KVector v1 = new KVector(5, 50);
        KVector v2 = new KVector();
        v2.parse("(5,50)");
        v2.parse("{5,50}");
        assertTrue(v1.equals(v2));
        v2.parse("[5,50]");
        assertTrue(v1.equals(v2));
        v2.parse("{(5,50)}");
        assertTrue(v1.equals(v2));
        v2.parse("[(5,50)]");
        assertTrue(v1.equals(v2));
        v2.parse("[{5,50}]");
        assertTrue(v1.equals(v2));
    }

    /**
     * Test applyBounds from KVector class.
     */
    @Test
    public void testApplyBounds() {
        /*
         * test if vt.x > lowx and vt.y > lowy (the result must be the same vt)
         */
        // init KVector v
        KVector v = new KVector(30, 30);
        // init KVector vt = v
        KVector vt = v;
        // init lower bound KVector v_lower_bound
        KVector v_lower_bound = new KVector(10, 10);
        // init upper bound KVector v_upper_bound
        KVector v_upper_bound = new KVector(40, 40);
        // test if vt.x > lowx (the result must be the same vt)
        vt.bound(v_lower_bound.x, v_lower_bound.y, v_upper_bound.x, v_upper_bound.y);
        assertTrue(vt.equals(v));
        
        /*
         * test if vt.x < lowx and vt.y < lowy (the result must be the vt(lowx,lowy))
         */
        // reinitialize vt
        vt = v;
        // init lower bound KVector v_lower_bound
        v_lower_bound = new KVector(40, 40);
        // init upper bound KVector v_upper_bound
        v_upper_bound = new KVector(60, 60);
        // test if vt.x < lowx and vt.y < lowy (the result must be the vt(lowx,lowy))
        vt.bound(v_lower_bound.x, v_lower_bound.y, v_upper_bound.x, v_upper_bound.y);
        assertTrue(vt.equals(v_lower_bound));
        
        /*
         * test if vt.x > highx and vt.y > highy (the result must be the vt(highx,highy))
         */
        // reinitialize vt
        vt = v;
        // init lower bound KVector v_lower_bound
        v_lower_bound = new KVector(20, 20);
        // init upper bound KVector v_upper_bound
        v_upper_bound = new KVector(30, 30);
        // test if vt.x > highx and vt.y > highy (the result must be the vt(highx,highy))
        vt.bound(v_lower_bound.x, v_lower_bound.y, v_upper_bound.x, v_upper_bound.y);
        assertTrue(vt.equals(v_upper_bound));
    }

}
