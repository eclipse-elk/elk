/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.spore.test;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import org.eclipse.elk.alg.common.compaction.oned.CompareFuzzy;
import org.eclipse.elk.alg.common.spore.ScanlineOverlapCheck;
import org.eclipse.elk.alg.spore.options.SporeMetaDataProvider;
import org.eclipse.elk.alg.spore.options.SporeOverlapRemovalOptions;
import org.eclipse.elk.core.RecursiveGraphLayoutEngine;
import org.eclipse.elk.core.data.LayoutMetaDataService;
import org.eclipse.elk.core.math.ElkMath;
import org.eclipse.elk.core.math.ElkRectangle;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.util.BasicProgressMonitor;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.eclipse.emf.ecore.util.EcoreUtil.Copier;
import org.junit.Test;

/**
 * Unit test for overlap removal using a scanline.
 */
public class ScanlineOverlapRemovalTest {
	
	/**
	 * Test whether the {@link ScanlineOverlapCheck} supports the overlap removal to find all
	 * overlaps in a special case.
	 */
	@Test
	public void scanlineTest() {
		// build test graphs
		ElkNode graph1 = ElkGraphUtil.createGraph();
		ElkNode n0 = ElkGraphUtil.createNode(graph1);
		ElkNode n1 = ElkGraphUtil.createNode(graph1);
		ElkNode n2 = ElkGraphUtil.createNode(graph1);
		ElkNode n3 = ElkGraphUtil.createNode(graph1);
		n0.setDimensions(160, 20);
		n1.setDimensions(160, 20);
		n2.setDimensions(20, 20);
		n3.setDimensions(20, 20);
		n0.setLocation(0, 30);
		n1.setLocation(150, 40);
		n2.setLocation(150, 0);
		n3.setLocation(150, 70);
		graph1.setProperty(CoreOptions.ALGORITHM, SporeOverlapRemovalOptions.ALGORITHM_ID);
		Copier copier = new Copier();
        ElkNode graph2 = (ElkNode) copier.copy(graph1);
        graph2.setProperty(SporeOverlapRemovalOptions.OVERLAP_REMOVAL_RUN_SCANLINE, false);
        
        // execute overlap removal with and without ScanlineOverlapCheck
        LayoutMetaDataService lService = LayoutMetaDataService.getInstance();
        lService.registerLayoutMetaDataProviders(new SporeMetaDataProvider());
        RecursiveGraphLayoutEngine lEngine = new RecursiveGraphLayoutEngine();
        lEngine.layout(graph1, new BasicProgressMonitor());
        lEngine.layout(graph2, new BasicProgressMonitor());
        
        // test
        assertFalse(hasOverlaps(graph1));
        assertTrue(hasOverlaps(graph2));
	}
	
	private boolean hasOverlaps(ElkNode graph) {
		
		for (ElkNode na : graph.getChildren()) {
			for (ElkNode nb : graph.getChildren()) {
				if (na != nb) {
					ElkRectangle r1 = new ElkRectangle(na.getX(), na.getY(), na.getWidth(), na.getHeight());
					ElkRectangle r2 = new ElkRectangle(nb.getX(), nb.getY(), nb.getWidth(), nb.getHeight());
					if (CompareFuzzy.lt(ElkMath.shortestDistance(r1, r2), 0)) {
						return true;
					}
				}
			}
		}
		
		return false;
	}
}
