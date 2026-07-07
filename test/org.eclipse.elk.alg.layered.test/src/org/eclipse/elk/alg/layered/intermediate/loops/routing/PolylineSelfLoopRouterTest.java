/*******************************************************************************
 * Copyright (c) 2019 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate.loops.routing;

import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.math.KVectorChain;
import org.junit.Assert;
import org.junit.Test;

/**
 * Tests whether the {@link PolylineSelfLoopRouter}'s corner cutting works as expected.
 */
public class PolylineSelfLoopRouterTest {
    
    private PolylineSelfLoopRouter router = new PolylineSelfLoopRouter();
    
    @Test
    public void testUsualCase() {
        KVectorChain input = new KVectorChain();
        input.add(new KVector(0, 0));
        input.add(new KVector(100, 0));
        input.add(new KVector(100, 100));
        input.add(new KVector(-100, 100));
        input.add(new KVector(-100, -100));
        input.add(new KVector(0, -100));
        
        KVectorChain expected = new KVectorChain();
        expected.add(new KVector(90, 0));
        expected.add(new KVector(100, 10));
        expected.add(new KVector(100, 90));
        expected.add(new KVector(90, 100));
        expected.add(new KVector(-90, 100));
        expected.add(new KVector(-100, 90));
        expected.add(new KVector(-100, -90));
        expected.add(new KVector(-90, -100));
        
        Assert.assertEquals(expected, router.cutCorners(input, 10));
    }
    
    @Test
    public void testSmallerDistance() {
        KVectorChain input = new KVectorChain();
        input.add(new KVector(0, 0));
        input.add(new KVector(100, 0));
        input.add(new KVector(100, 100));
        input.add(new KVector(-100, 100));
        input.add(new KVector(-100, -100));
        input.add(new KVector(0, -100));
        
        KVectorChain expected = new KVectorChain();
        expected.add(new KVector(50, 0));
        expected.add(new KVector(100, 50));
        expected.add(new KVector(100, 50));
        expected.add(new KVector(50, 100));
        expected.add(new KVector(-20, 100));
        expected.add(new KVector(-100, 20));
        expected.add(new KVector(-100, -50));
        expected.add(new KVector(-50, -100));
        
        Assert.assertEquals(expected, router.cutCorners(input, 80));
    }
    
}
