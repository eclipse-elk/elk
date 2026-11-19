/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.

 * This Source Code may also be made available under the following Secondary
 * Licenses when the conditions for such availability set forth in the Eclipse
 * Public License v. 2.0 are satisfied: GPL-3.0 which is available at
 * https://www.gnu.org/licenses/gpl-3.0-standalone.html.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.spore.test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Arrays;
import java.util.Map;
import java.util.Set;

import org.eclipse.elk.alg.common.compaction.oned.CompareFuzzy;
import org.eclipse.elk.alg.common.spore.Node;
import org.eclipse.elk.alg.common.spore.ScanlineOverlapCheck;
import org.eclipse.elk.alg.common.utils.SVGImage;
import org.eclipse.elk.alg.spore.options.SporeMetaDataProvider;
import org.eclipse.elk.alg.spore.options.SporeOverlapRemovalOptions;
import org.eclipse.elk.core.RecursiveGraphLayoutEngine;
import org.eclipse.elk.core.data.LayoutMetaDataService;
import org.eclipse.elk.core.math.ElkMath;
import org.eclipse.elk.core.math.ElkRectangle;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.util.BasicProgressMonitor;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.eclipse.emf.ecore.util.EcoreUtil.Copier;
import org.junit.Test;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

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
	
	/**
	 * Overlap removal must not fail with an {@link IllegalStateException} from the
	 * {@link ScanlineOverlapCheck} on three overlapping nodes, and must remove their overlaps.
	 */
	@Test
	public void scanlineThreeOverlappingNodesTest() {
		ElkNode graph = ElkGraphUtil.createGraph();
		ElkNode n0 = ElkGraphUtil.createNode(graph);
		ElkNode n1 = ElkGraphUtil.createNode(graph);
		ElkNode n2 = ElkGraphUtil.createNode(graph);
		n0.setDimensions(40, 30);
		n1.setDimensions(10, 20);
		n2.setDimensions(10, 30);
		n0.setLocation(10, 10);
		n1.setLocation(10, 0);
		n2.setLocation(20, 5);
		graph.setProperty(CoreOptions.ALGORITHM, SporeOverlapRemovalOptions.ALGORITHM_ID);
		assertTrue(hasOverlaps(graph));

		LayoutMetaDataService.getInstance().registerLayoutMetaDataProviders(new SporeMetaDataProvider());
		new RecursiveGraphLayoutEngine().layout(graph, new BasicProgressMonitor());

		assertFalse(hasOverlaps(graph));
	}
	
	/**
	 * The {@link ScanlineOverlapCheck} must report an overlap with a node that lies, in the order of
	 * left sides, behind a node that does not overlap the inserted node: here the narrow node b
	 * between the wide node a and the node c, which both overlap the inserted node p.
	 */
	@Test
	public void scanlineOverlapBehindNarrowNodeTest() {
		Map<String, Node> nodes = Maps.newLinkedHashMap();
		nodes.put("a", node(0, 0, 100, 10));
		nodes.put("b", node(10, 0, 5, 10));
		nodes.put("c", node(55, 0, 10, 10));
		nodes.put("p", node(50, 5, 10, 10));
		
		assertEquals(Sets.newTreeSet(Arrays.asList("a-b", "a-c", "a-p", "c-p")), findOverlaps(nodes));
	}
	
	/**
	 * The {@link ScanlineOverlapCheck} must report an overlap with a node whose left side lies right of
	 * the inserted node's left side, after an overlap with a node left of the inserted node.
	 */
	@Test
	public void scanlineOverlapRightOfInsertedNodeTest() {
		Map<String, Node> nodes = Maps.newLinkedHashMap();
		nodes.put("a", node(0, 0, 100, 10));
		nodes.put("c", node(55, 0, 10, 10));
		nodes.put("p", node(50, 5, 10, 10));
		
		assertEquals(Sets.newTreeSet(Arrays.asList("a-c", "a-p", "c-p")), findOverlaps(nodes));
	}
	
	/**
	 * Sweeps the given nodes and returns the overlapping pairs the {@link ScanlineOverlapCheck} reports,
	 * each as the two node names in alphabetical order, joined by a dash.
	 */
	private Set<String> findOverlaps(Map<String, Node> nodes) {
		Map<Node, String> names = Maps.newHashMap();
		nodes.forEach((name, node) -> names.put(node, name));
		Set<String> overlaps = Sets.newTreeSet();
		new ScanlineOverlapCheck((n1, n2) -> {
			String name1 = names.get(n1);
			String name2 = names.get(n2);
			overlaps.add(name1.compareTo(name2) < 0 ? name1 + "-" + name2 : name2 + "-" + name1);
		}, new SVGImage(null)).sweep(Lists.newArrayList(nodes.values()));
		return overlaps;
	}
	
	private Node node(double x, double y, double width, double height) {
		return new Node(new KVector(x + width / 2, y + height / 2), new ElkRectangle(x, y, width, height));
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
