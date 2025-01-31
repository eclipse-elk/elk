package org.eclipse.elk.core.math;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

public class KVectorChainTest {

    /**
     * Test parse from KVectorChain class.
     */
    @Test
    public void testParse() {
        KVector v0 = new KVector(5, 50);
        KVector v1 = new KVector(10, 50);
        KVector v2 = new KVector(30, 50);
        KVectorChain kv = new KVectorChain();
        kv.parse("{(5,50),(10,50),(30,50)}");

        assertTrue(v0.equals(kv.get(0)));
        assertTrue(v1.equals(kv.get(1)));
        assertTrue(v2.equals(kv.get(2)));
        
        // ignore a dangling element
        kv = new KVectorChain();
        kv.parse("{(5,50),(10,50),(30,)}");

        assertTrue(v0.equals(kv.get(0)));
        assertTrue(v1.equals(kv.get(1)));
        assertTrue(kv.size() == 2);
        
        // some weird syntax
        kv = new KVectorChain();
        kv.parse("{(5; 50 ], [10 , 50 ),(30,,,)}");

        assertTrue(v0.equals(kv.get(0)));
        assertTrue(v1.equals(kv.get(1)));
        assertTrue(kv.size() == 2);
    }

    /** 
     * Tests the IllegalArgumentException of parse from KVectorChain class.
     */
    @Test(expected = IllegalArgumentException.class)
    public void testParseIllegalArgumentException() {
        KVectorChain kv = new KVectorChain();
        kv.parse("{(5,a),(10,50),(30,50)}");
    }

    /**
     * Test getLength from KVectorChain class. This function test if the length result ist 0 when
     * KVectorChain composed of 3 overlapping KVectors. And test for 3 differing KVectors.
     */
    @Test
    public void testGetLength() {
        // 3 overlaped KVectors
        KVectorChain kv = new KVectorChain();
        kv.parse("{(10,50),(10,50),(10,50)}");
        assertEquals(0, kv.totalLength(), 0);
        
        // 3 differing KVecors
        kv.parse("{(10,0),(10,20),(10,30)}");
        assertEquals(30, kv.totalLength(), 0);
    }

    /**
     * Test getPointOnLine from KVectorChain class.
     */
    @Test
    public void testGetPointOnLine() {
        KVector v0 = new KVector(5, 50);
        KVector v1 = new KVector(10, 50);
        KVector v2 = new KVector(30, 50);
        KVectorChain kv = new KVectorChain(v0, v1, v2);
        
        // test if resturns v0 for distance = 0
        assertTrue(v0.equals(kv.pointOnLine(0)));
        
        // test if resturns v1 for distance = 5
        assertTrue(v1.equals(kv.pointOnLine(5)));
        
        // test if resturns endpoint  for distance > KVectorChain's length
        assertTrue(v2.equals(kv.pointOnLine(40)));
    }

}
