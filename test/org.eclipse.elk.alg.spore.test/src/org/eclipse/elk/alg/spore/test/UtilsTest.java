/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.spore.test;

import static org.junit.Assert.assertTrue;

import java.util.List;

import org.eclipse.elk.alg.common.compaction.oned.CompareFuzzy;
import org.eclipse.elk.alg.common.spore.Node;
import org.eclipse.elk.alg.common.utils.Utils;
import org.eclipse.elk.core.math.ElkMath;
import org.eclipse.elk.core.math.ElkRectangle;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.util.Pair;
import org.junit.Before;
import org.junit.Test;

import com.google.common.collect.Lists;

/**
 * Unit test for spore utils.
 */
public class UtilsTest {

	/**
	 * Testing the computed overlap of two rectangles by moving one accordingly and checking 
	 * whether the distance becomes >= 0.
	 */
	@Test
	public void overlapTest() {
		ElkRectangle r1 = new ElkRectangle(0, 0, 40, 80);
        ElkRectangle r2 = new ElkRectangle(0, 0, 70, 20);
        ElkRectangle r3 = new ElkRectangle(10, 50, 70, 20);
        ElkRectangle r4 = new ElkRectangle(-40, 30, 70, 20);
        ElkRectangle r5 = new ElkRectangle(-60, 70, 70, 20);
        ElkRectangle r6 = new ElkRectangle(-10, 70, 70, 20);
        ElkRectangle r7 = new ElkRectangle(-20, -10, 70, 20);
        ElkRectangle r8 = new ElkRectangle(10, 20, 20, 20);
        ElkRectangle r9 = new ElkRectangle(-20, -20, 100, 120);
        ElkRectangle r10 = new ElkRectangle(0, -0.001, 40, 80);

        List<ElkRectangle> rs = Lists.newArrayList(r2,r3,r4,r5,r6,r7,r8,r9,r10);
        rs.forEach(r -> assertTrue(testOverlapComputation(r1, r)));
	}
	
	private boolean testOverlapComputation(ElkRectangle r1, ElkRectangle r2) {
		// compute position for r2 to not overlap r1
		double o = Utils.overlap(r1, r2);
        KVector c1 = r1.getCenter();
        KVector c2 = r2.getCenter();
        KVector d = c2.clone().sub(c1);
        KVector c3 = c1.clone().add(d.scale(o));
        ElkRectangle r3 = new ElkRectangle(r2);
        r3.move(c3.clone().sub(c2));
        
        return CompareFuzzy.ge(ElkMath.shortestDistance(r1, r3), 0.0);
	}
	
	private List<ElkRectangle> rectangles;
	private ElkRectangle r1;
	@Before
	public void init() {
		r1 = new ElkRectangle(0, 0, 20, 60);
		ElkRectangle r2 = new ElkRectangle(40, 20, 20, 20);
		ElkRectangle r3 = new ElkRectangle(40, 40, 20, 20);
		ElkRectangle r4 = new ElkRectangle(30, 70, 20, 20);
		ElkRectangle r5 = new ElkRectangle(20, 80, 20, 20);
		ElkRectangle r6 = new ElkRectangle(10, 80, 20, 20);
		ElkRectangle r7 = new ElkRectangle(0, 80, 20, 20);
		ElkRectangle r8 = new ElkRectangle(-30, 70, 20, 20);
		ElkRectangle r9 = new ElkRectangle(-40, 40, 20, 20);
		ElkRectangle r10 = new ElkRectangle(-40, 20, 20, 20);
		ElkRectangle r11 = new ElkRectangle(-30, -20, 20, 20);
		ElkRectangle r12 = new ElkRectangle(-20, -40, 20, 20);
		ElkRectangle r13 = new ElkRectangle(0, -30, 20, 20);
		ElkRectangle r14 = new ElkRectangle(30, -30, 20, 20);
		ElkRectangle r15 = new ElkRectangle(20, 0, 20, 20);
		ElkRectangle r16 = new ElkRectangle(0, 60, 20, 20);
		
		rectangles = Lists.newArrayList(r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13,r14,r15,r16);
	}
	
	/**
	 * Testing the computed underlap and distance of two nodes by moving one accordingly and checking 
	 * whether the distance becomes 0.
	 */
	@Test
	public void underlapTest() {
		rectangles.forEach(r -> assertTrue(testUnderlapComputation(r1, r)));
		
	}
	
	private boolean testUnderlapComputation(ElkRectangle r1, ElkRectangle r2) {
		Node n1 = new Node(r1.getCenter(), r1);
		Node n2 = new Node(r2.getCenter(), r2);
		
		double underlap = n1.underlap(n2);
		// Additionally, the underlap should be equal to the distance in the direction
		// given by the nodes' centers.
		assertTrue(CompareFuzzy.eq(underlap, n1.distance(n2, n1.vertex.clone().sub(n2.vertex))));
		
		// move + check shortest distance
		n2.translate(n1.vertex.clone().sub(n2.vertex).scaleToLength(underlap));
		
		return CompareFuzzy.eq(ElkMath.shortestDistance(n1.rect, n2.rect), 0.0);
	}
	
	/**
	 * Distance test using different directions.
	 */
	@Test
	public void distanceTest() {
		List<Pair<KVector, Boolean>> vectors = Lists.newArrayList();
		vectors.add(Pair.of(new KVector(-20, 20), true)); // true if the nodes should collide
		vectors.add(Pair.of(new KVector(-80, 0), false));
		vectors.add(Pair.of(new KVector(-20, 9), false));
		vectors.add(Pair.of(new KVector(0, 50), false));
		vectors.add(Pair.of(new KVector(-9.99, 50), false));
		vectors.add(Pair.of(new KVector(60, 60), false));
		vectors.add(Pair.of(new KVector(-30, 50), true));
		vectors.add(Pair.of(new KVector(-20, 130), true));
		vectors.add(Pair.of(new KVector(-20, -21), false));
		
//		Arrays.stream(new File(ElkUtil.debugFolderPath("test")).listFiles()).forEach(File::delete);
		
		Node n1 = new Node(r1.getCenter(), r1);
		for (Pair<KVector, Boolean> v : vectors) {
			Node n2 = new Node(rectangles.get(12).getCenter(), new ElkRectangle(rectangles.get(12)));
			
			double distance = n1.distance(n2, v.getFirst());
			
			// if they should collide, move them, otherwise check whether the returned distance is infinite
			if (v.getSecond()) {
				// move + check shortest distance
				n2.translate(v.getFirst().scaleToLength(distance));
				assertTrue(CompareFuzzy.eq(ElkMath.shortestDistance(n1.rect, n2.rect), 0.0));
				
			} else {
				assertTrue(Double.isInfinite(distance));
			}
		}
	}
}
