/*******************************************************************************
 * Copyright (c) 2017, 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.core.math;

import static org.junit.Assert.*;
import org.junit.Test;

/**
 * Various tests for the {@link ElkMath} class.
 */
public class ElkMathTest {

    @Test
    public void testPointContains() {
        ElkRectangle rect = new ElkRectangle(23, 14, 20, 20);
        KVector contained = new KVector(26, 19);
        assertTrue(ElkMath.contains(rect, contained));
        
        KVector notContained = new KVector(10, 9); 
        assertTrue(!ElkMath.contains(rect, notContained));        
                
        KVector onBorder = new KVector(23, 20);
        assertTrue(!ElkMath.contains(rect, onBorder));
        
        KVector onCorner = new KVector(23, 14);
        assertTrue(!ElkMath.contains(rect, onCorner));
    }
    
    @Test 
    public void testLineContains() {
        ElkRectangle rect = new ElkRectangle(23, 14, 20, 20);
        KVector line11 = new KVector(24, 20);
        KVector line12 = new KVector(40, 32);
        assertTrue(ElkMath.contains(rect, line11, line12));
        
        KVector line21 = new KVector(10, 10); // outside
        KVector line22 = new KVector(40, 32);
        assertTrue(!ElkMath.contains(rect, line21, line22));
    }
    
    @Test
    public void testPathContains() {
        ElkRectangle rect = new ElkRectangle(23, 14, 20, 20);
        
        KVectorChain path = new KVectorChain();
        path.add(24, 15);
        path.add(27, 20);
        path.add(39, 30);
        path.add(29, 19);
        assertTrue(ElkMath.contains(rect, path));
        
        // on border
        KVectorChain path2 = new KVectorChain(path);
        path2.add(23, 14);
        assertTrue(!ElkMath.contains(rect, path2));
        
        // outside
        path.add(10, 10);
        assertTrue(!ElkMath.contains(rect, path));
    }
    
    @Test
    public void testLineLineIntersect() {
        KVector l11 = new KVector(10, 10);
        KVector l12 = new KVector(20, 20);
        
        // cross
        KVector l21 = new KVector(11, 21);
        KVector l22 = new KVector(21, 11);
        assertTrue(ElkMath.intersects(l11, l12, l21, l22));
        
        // touch
        l21 = new KVector(10, 10);
        l22 = new KVector(15, 10);
        assertTrue(!ElkMath.intersects(l11, l12, l21, l22));
        
        // no cross
        l21 = new KVector(1, 2);
        l22 = new KVector(2, 1);
        assertTrue(!ElkMath.intersects(l11, l12, l21, l22));
        
        // same line
        assertTrue(!ElkMath.intersects(l11, l12, l11.clone(), l12.clone()));
        assertTrue(!ElkMath.intersects(l11, l12, l12.clone(), l11.clone()));
        
        // parallel
        l21 = new KVector(11, 1);
        l22 = new KVector(21, 21);
        assertTrue(!ElkMath.intersects(l11, l12, l21, l22));
    }
    
    @Test
    public void testPathIntersects() {
        ElkRectangle rect = new ElkRectangle(23, 14, 20, 20);
        
        KVectorChain path = new KVectorChain();
        path.add(24, 15);
        path.add(27, 20);
        path.add(39, 30);
        path.add(29, 19);
        assertTrue(!ElkMath.intersects(rect, path));
        
        // on border
        KVectorChain path2 = new KVectorChain(path);
        path2.add(23, 14);
        assertTrue(!ElkMath.intersects(rect, path2));
        
        // cross
        path.add(10, 10);
        assertTrue(ElkMath.intersects(rect, path));
    }

    /**
     * Tests some valid combinations of Factl from ElkMath class.
     */
    @Test
    public void testFactl() {
        assertEquals(1, ElkMath.factl(0));
        assertEquals(1, ElkMath.factl(1));
        assertEquals(2432902008176640000L, ElkMath.factl(20));
    }

    /**
     * Tests the Little argument exception of Factl from ElkMath class.
     */
    @Test(expected = IllegalArgumentException.class)
    public void testFactlLittleIllegalArgumentException() {
        ElkMath.factl(-50);
    }

    /**
     * Tests the Big argument exception of Factl from ElkMath class.
     */
    @Test(expected = IllegalArgumentException.class)
    public void testFactlBigIllegalArgumentException() {
        ElkMath.factl(21);
    }

    /**
     * Tests some valid combinations of Factd from ElkMath class.
     */
    @Test
    public void testFactd() {
        assertEquals(1, ElkMath.factd(0), 0.0);
        assertEquals(1, ElkMath.factd(1), 0.0);
    }

    /**
     * Tests the Little argument exception of Factd from ElkMath class.
     */
    @Test(expected = IllegalArgumentException.class)
    public void testFacdlLittleIllegalArgumentException() {
        ElkMath.factd(-1);
    }

    /**
     * Tests some valid combinations of binomiall from ElkMath class.
     */
    @Test
    public void testBinomiall() {
        assertEquals(1, ElkMath.binomiall(2, 0));
        assertEquals(1, ElkMath.binomiall(20, 20));
        assertEquals(2, ElkMath.binomiall(2, 1));
    }

    /**
     * Tests the Little argument exception of binomiall from ElkMath class.
     */
    @Test(expected = IllegalArgumentException.class)
    public void testBinomialllLittleIllegalArgumentException() {
        ElkMath.binomiall(-1, 1);
    }

    /**
     * Tests some valid combinations of binomiald from ElkMath class.
     */
    @Test
    public void testBinomiald() {
        assertEquals(1, ElkMath.binomiald(2, 0), 0);
        assertEquals(1, ElkMath.binomiald(20, 20), 0);
        assertEquals(2, ElkMath.binomiald(2, 1), 0);
    }

    /**
     * Tests the Little argument exception of binomiald from ElkMath class.
     */
    @Test(expected = IllegalArgumentException.class)
    public void testBinomialdLittleIllegalArgumentException() {
        ElkMath.binomiald(-1, 1);
    }

    /**
     * Tests some valid combinations of pow from ElkMath class.
     */
    @Test
    public void testPow() {
        double ad = 10;
        float af = 10;
        assertEquals(1, ElkMath.powd(ad, 0), 0);
        assertEquals(1, ElkMath.powf(af, 0), 0);
        assertEquals(100, ElkMath.powd(ad, 2), 0);
        assertEquals(100, ElkMath.powf(af, 2), 0);
    }

    /**
     * Tests some valid combinations of CalcBezierPoints from ElkMath class.
     */
    @Test
    public void testCalcBezierPoints() {
        // some KVectors
        KVector kvector1 = new KVector(10, 10);
        KVector kvector2 = new KVector(20, 20);
        KVector kvector3 = new KVector(30, 30);
        KVector kvector4 = new KVector(50, 50);
        
        // test if the last KVector of the result similar to kvector4
        KVector[] result = ElkMath.approximateBezierSegment(20, kvector1, kvector2, kvector3, kvector4);
        assertEquals(kvector4.x, result[result.length - 1].x, 0.000000001);
        assertEquals(kvector4.y, result[result.length - 1].y, 0.000000001);

        // some KVectors with y=10
        kvector1 = new KVector(50, 10);
        kvector2 = new KVector(70, 10);
        kvector3 = new KVector(80, 10);
        kvector4 = new KVector(100, 10);
        
        // test if the all result-KVectors with y=10
        result = ElkMath.approximateBezierSegment(20, kvector1, kvector2, kvector3, kvector4);
        for (KVector k : result) {
            assertEquals(10, k.y, 0.000000001);
        }
    }

    /**
     * Tests some valid combinations of appoximateSpline from ElkMath class.
     */
    @Test
    public void testAppoximateSpline() {
        // some KVectors
        KVector kvector1 = new KVector(10, 10);
        KVector kvector2 = new KVector(20, 20);
        KVector kvector3 = new KVector(30, 30);
        KVector kvector4 = new KVector(50, 50);
        
        // test if the last KVector of the result similar to kvector4
        KVector[] vectors = ElkMath.approximateBezierSegment(20, kvector1, kvector2, kvector3, kvector4);
        KVectorChain controlPoints = new KVectorChain(vectors);
        KVectorChain result = ElkMath.approximateBezierSpline(controlPoints);
        KVector k = result.get(result.size() - 1);
        assertEquals(kvector4.x, k.x, 0.000000001);
        assertEquals(kvector4.y, k.y, 0.000000001);

        // some KVectors with y=10
        kvector1 = new KVector(50, 10);
        kvector2 = new KVector(70, 10);
        kvector3 = new KVector(80, 10);
        kvector4 = new KVector(100, 10);
        
        // test if the all result-KVectors with y=10
        vectors = ElkMath.approximateBezierSegment(20, kvector1, kvector2, kvector3, kvector4);
        controlPoints = new KVectorChain(vectors);
        result = ElkMath.approximateBezierSpline(controlPoints);

        for (KVector kv : result) {
            assertEquals(10, kv.y, 0.000000001);
        }
    }

    /**
     * Tests some valid combinations of distanceFromSpline from ElkMath class.
     */
    @Test
    public void testDistanceFromSpline() {
        // some KVectors
        KVector kvector1 = new KVector(10, 10);
        KVector kvector2 = new KVector(20, 20);
        KVector kvector3 = new KVector(30, 30);
        KVector kvector4 = new KVector(50, 50);
        
        // test if the result is 0 when kvector4 = needle
        KVector needle = kvector4;
        double result = ElkMath.distanceFromBezierSegment(kvector1, kvector2, kvector3, kvector4,
                needle);
        assertEquals(0, result, 0.01);
        
        // test if the result is 0 when kvector3 = needle
        needle = kvector3;
        result = ElkMath.distanceFromBezierSegment(kvector1, kvector2, kvector3, kvector4, needle);
        assertEquals(0, result, 0.01);
        
        // test if the result is 0 when kvector2 = needle
        needle = kvector2;
        result = ElkMath.distanceFromBezierSegment(kvector1, kvector2, kvector3, kvector4, needle);
        assertEquals(0, result, 0.01);
        
        // test if the result is 0 when kvector1 = needle
        needle = kvector1;
        result = ElkMath.distanceFromBezierSegment(kvector1, kvector2, kvector3, kvector4, needle);
        assertEquals(0, result, 0.01);
    }

    /**
     * Tests some valid combinations of maxi from ElkMath class.
     */
    @Test
    public void testMax() {
        // test if the max is 7
        assertEquals(7, ElkMath.maxi(1, 7, 5, 6));
        assertEquals(7, ElkMath.maxf(1, 7, 5, 6), 0);
        assertEquals(7, ElkMath.maxd(1, 7, 5, 6), 0);
    }

    /**
     * Tests some valid combinations of mini,minf,mind from ElkMath class.
     */
    @Test
    public void testMin() {
        // test if the mini is 1
        assertEquals(1, ElkMath.mini(1, 7, 5, 6));
        // test if the mini is 0
        assertEquals(0, ElkMath.mini(8, 1, 9, 0));
        // test if the mini is 8
        assertEquals(8, ElkMath.mini(8, 8, 8, 8));
        // test if the minf is 1
        assertEquals(1, ElkMath.minf(1, 7, 5, 6), 0);
        // test if the minf is 0
        assertEquals(0, ElkMath.minf(8, 1, 9, 0), 0);
        // test if the minf is 8
        assertEquals(8, ElkMath.minf(8, 8, 8, 8), 0);
        // test if the mind is 1
        assertEquals(1, ElkMath.mind(1, 7, 5, 6), 0);
        // test if the mind is 0
        assertEquals(0, ElkMath.mind(8, 1, 9, 0), 0);
        // test if the mind is 8
        assertEquals(8, ElkMath.mind(8, 8, 8, 8), 0);
    }

    /**
     * Tests some valid combinations of averagei from ElkMath class.
     */
    @Test
    public void testAverage() {
        // test if the averagei is 4
        assertEquals(4, ElkMath.averagel((long) 5, (long) 8, (long) 2, (long) 1));
        // test if the averagei is 2
        assertEquals(2, ElkMath.averagel((long) 5, (long) 0, (long) 2, (long) 1));
        // test if the averagef is 4
        assertEquals(4, ElkMath.averagef(5, 8, 2, 1), 0);
        // test if the averagef is 2
        assertEquals(2, ElkMath.averagef(5, 0, 2, 1), 0);
        // test if the averaged is 4
        assertEquals(4, ElkMath.averaged(5, 8, 2, 1), 0);
        // test if the averaged is 2
        assertEquals(2, ElkMath.averaged(5, 0, 2, 1), 0);
    }
    
}
