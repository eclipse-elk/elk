/*******************************************************************************
 * Copyright (c) 2026 Kiel University and others.
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
import static org.junit.Assert.assertNotEquals;

import java.util.ArrayList;
import java.util.List;

import org.eclipse.elk.alg.spore.options.SporeCompactionOptions;
import org.eclipse.elk.alg.spore.options.SporeMetaDataProvider;
import org.eclipse.elk.alg.spore.options.SporeOverlapRemovalOptions;
import org.eclipse.elk.core.RecursiveGraphLayoutEngine;
import org.eclipse.elk.core.data.LayoutMetaDataService;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.util.BasicProgressMonitor;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.junit.BeforeClass;
import org.junit.Test;

/**
 * Tests that the SPOrE algorithms shift nodes with identical centers apart reproducibly, controlled by the
 * {@link CoreOptions#RANDOM_SEED random seed} option.
 */
public class RandomSeedTest {

	/** How often each graph is laid out to check that the result does not change. */
	private static final int RUNS = 5;
	/** Number of nodes that share one center point. */
	private static final int COINCIDENT_NODES = 3;
	/** Width and height of each node. */
	private static final double NODE_SIZE = 30;

	@BeforeClass
	public static void registerSpore() {
		LayoutMetaDataService.getInstance().registerLayoutMetaDataProviders(new SporeMetaDataProvider());
	}

	@Test
	public void overlapRemovalIsReproducibleWithDefaultSeed() {
		assertReproducible(SporeOverlapRemovalOptions.ALGORITHM_ID, null);
	}

	@Test
	public void overlapRemovalIsReproducibleWithExplicitSeed() {
		assertReproducible(SporeOverlapRemovalOptions.ALGORITHM_ID, 42);
	}

	@Test
	public void compactionIsReproducibleWithDefaultSeed() {
		assertReproducible(SporeCompactionOptions.ALGORITHM_ID, null);
	}

	@Test
	public void compactionIsReproducibleWithExplicitSeed() {
		assertReproducible(SporeCompactionOptions.ALGORITHM_ID, 42);
	}

	@Test
	public void defaultSeedIsOne() {
		assertEquals(layoutCoincidentNodes(SporeOverlapRemovalOptions.ALGORITHM_ID, 1),
				layoutCoincidentNodes(SporeOverlapRemovalOptions.ALGORITHM_ID, null));
	}

	@Test
	public void seedSelectsTheShift() {
		assertNotEquals(layoutCoincidentNodes(SporeOverlapRemovalOptions.ALGORITHM_ID, 1),
				layoutCoincidentNodes(SporeOverlapRemovalOptions.ALGORITHM_ID, 2));
	}

	/**
	 * Lays out the same graph several times and asserts that every run places the nodes at the same positions.
	 */
	private void assertReproducible(final String algorithmId, final Integer seed) {
		List<KVector> first = layoutCoincidentNodes(algorithmId, seed);
		for (int run = 1; run < RUNS; run++) {
			assertEquals("run " + run, first, layoutCoincidentNodes(algorithmId, seed));
		}
	}

	/**
	 * Lays out a graph whose nodes all have the same center and returns their resulting positions.
	 *
	 * @param algorithmId
	 *            the SPOrE algorithm to apply.
	 * @param seed
	 *            the random seed to set, or {@code null} to leave the option at its default.
	 * @return the position of each node, in the order of the graph's children.
	 */
	private List<KVector> layoutCoincidentNodes(final String algorithmId, final Integer seed) {
		ElkNode graph = ElkGraphUtil.createGraph();
		graph.setProperty(CoreOptions.ALGORITHM, algorithmId);
		if (seed != null) {
			graph.setProperty(CoreOptions.RANDOM_SEED, seed);
		}
		for (int i = 0; i < COINCIDENT_NODES; i++) {
			ElkNode node = ElkGraphUtil.createNode(graph);
			node.setDimensions(NODE_SIZE, NODE_SIZE);
		}

		new RecursiveGraphLayoutEngine().layout(graph, new BasicProgressMonitor());

		List<KVector> positions = new ArrayList<>();
		for (ElkNode node : graph.getChildren()) {
			positions.add(new KVector(node.getX(), node.getY()));
		}
		return positions;
	}

}
