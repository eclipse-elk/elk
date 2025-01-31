/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.core.alg;

import static org.junit.Assert.assertEquals;

import org.junit.Ignore;
import org.junit.Test;

/**
 * Tests for {@link LayoutProcessorConfiguration}.
 */
public class LayoutProcessorConfigurationTest {

	/**
	 * Test method for {@link org.eclipse.elk.alg.common.structuring.LayoutProcessorConfiguration#clear()}.
	 */
	@Test
	@Ignore
	public void testClear() {
		LayoutProcessorConfiguration<TestPhases, StringBuffer> config = LayoutProcessorConfiguration.create();
		config.before(TestPhases.PHASE_1).add(TestProcessors.PROCESSOR_1);
		// FIXME this method is package-visible
//		config.clear();
		
		assertEquals(0, config.processorsBefore(TestPhases.PHASE_1).size());
		assertEquals(0, config.processorsAfter(TestPhases.PHASE_1).size());
		assertEquals(0, config.processorsBefore(TestPhases.PHASE_2).size());
		assertEquals(0, config.processorsAfter(TestPhases.PHASE_2).size());
	}

	/**
	 * Test method for {@link org.eclipse.elk.alg.common.structuring.LayoutProcessorConfiguration#before(java.lang.Enum)}.
	 */
	@Test
	public void testBefore() {
		LayoutProcessorConfiguration<TestPhases, StringBuffer> config = LayoutProcessorConfiguration.create();
		config.before(TestPhases.PHASE_1)
			.add(TestProcessors.PROCESSOR_1)
		.before(TestPhases.PHASE_2)
			.add(TestProcessors.PROCESSOR_2)
			.add(TestProcessors.PROCESSOR_3);
		
		assertEquals(1, config.processorsBefore(TestPhases.PHASE_1).size());
		assertEquals(2, config.processorsAfter(TestPhases.PHASE_1).size());
		assertEquals(2, config.processorsBefore(TestPhases.PHASE_2).size());
		assertEquals(0, config.processorsAfter(TestPhases.PHASE_2).size());
	}

	/**
	 * Test method for {@link org.eclipse.elk.alg.common.structuring.LayoutProcessorConfiguration#after(java.lang.Enum)}.
	 */
	@Test
	public void testAfter() {
		LayoutProcessorConfiguration<TestPhases, StringBuffer> config = LayoutProcessorConfiguration.create();
		config.after(TestPhases.PHASE_1)
			.add(TestProcessors.PROCESSOR_1)
		.after(TestPhases.PHASE_2)
			.add(TestProcessors.PROCESSOR_2)
			.add(TestProcessors.PROCESSOR_3);
		
		assertEquals(0, config.processorsBefore(TestPhases.PHASE_1).size());
		assertEquals(1, config.processorsAfter(TestPhases.PHASE_1).size());
		assertEquals(1, config.processorsBefore(TestPhases.PHASE_2).size());
		assertEquals(2, config.processorsAfter(TestPhases.PHASE_2).size());
	}

	/**
	 * Test method for {@link org.eclipse.elk.alg.common.structuring.LayoutProcessorConfiguration#addBefore(java.lang.Enum, org.eclipse.elk.core.alg.ILayoutProcessorFactory)}.
	 */
	@Test
	public void testAddBefore1() {
		LayoutProcessorConfiguration<TestPhases, StringBuffer> config = LayoutProcessorConfiguration.create();
		config.before(TestPhases.PHASE_1);
		config.addBefore(TestPhases.PHASE_2, TestProcessors.PROCESSOR_1);
		
		assertEquals(0, config.processorsBefore(TestPhases.PHASE_1).size());
		assertEquals(1, config.processorsAfter(TestPhases.PHASE_1).size());
		assertEquals(1, config.processorsBefore(TestPhases.PHASE_2).size());
		assertEquals(0, config.processorsAfter(TestPhases.PHASE_2).size());
	}

    /**
     * Test method for {@link org.eclipse.elk.alg.common.structuring.LayoutProcessorConfiguration#addBefore(java.lang.Enum, org.eclipse.elk.core.alg.ILayoutProcessorFactory)}.
     */
    @Test(expected = IllegalStateException.class)
    public void testAddBefore2() {
        LayoutProcessorConfiguration<TestPhases, StringBuffer> config = LayoutProcessorConfiguration.create();
        config.before(TestPhases.PHASE_1);
        config.addBefore(TestPhases.PHASE_2, TestProcessors.PROCESSOR_1);
        config.add(TestProcessors.PROCESSOR_2);
    }

	/**
	 * Test method for {@link org.eclipse.elk.alg.common.structuring.LayoutProcessorConfiguration#addAfter(java.lang.Enum, org.eclipse.elk.core.alg.ILayoutProcessorFactory)}.
	 */
	@Test
	public void testAddAfter1() {
		LayoutProcessorConfiguration<TestPhases, StringBuffer> config = LayoutProcessorConfiguration.create();
		config.before(TestPhases.PHASE_1);
		config.addAfter(TestPhases.PHASE_1, TestProcessors.PROCESSOR_1);
		
		assertEquals(0, config.processorsBefore(TestPhases.PHASE_1).size());
		assertEquals(1, config.processorsAfter(TestPhases.PHASE_1).size());
		assertEquals(1, config.processorsBefore(TestPhases.PHASE_2).size());
		assertEquals(0, config.processorsAfter(TestPhases.PHASE_2).size());
	}

    /**
     * Test method for {@link org.eclipse.elk.alg.common.structuring.LayoutProcessorConfiguration#addAfter(java.lang.Enum, org.eclipse.elk.core.alg.ILayoutProcessorFactory)}.
     */
    @Test(expected = IllegalStateException.class)
    public void testAddAfter2() {
        LayoutProcessorConfiguration<TestPhases, StringBuffer> config = LayoutProcessorConfiguration.create();
        config.before(TestPhases.PHASE_1);
        config.addAfter(TestPhases.PHASE_1, TestProcessors.PROCESSOR_1);
        config.add(TestProcessors.PROCESSOR_2);
    }

	/**
	 * Test method for {@link org.eclipse.elk.alg.common.structuring.LayoutProcessorConfiguration#addAll(org.eclipse.elk.alg.common.structuring.LayoutProcessorConfiguration)}.
	 */
	@Test
	public void testAddAll() {
		LayoutProcessorConfiguration<TestPhases, StringBuffer> config = LayoutProcessorConfiguration.create();
		config.after(TestPhases.PHASE_1)
			.add(TestProcessors.PROCESSOR_1)
		.after(TestPhases.PHASE_2)
			.add(TestProcessors.PROCESSOR_2)
			.add(TestProcessors.PROCESSOR_3);
		
		LayoutProcessorConfiguration<TestPhases, StringBuffer> config2 = LayoutProcessorConfiguration.create();
		config2.addAll(config);

		assertEquals(0, config2.processorsBefore(TestPhases.PHASE_1).size());
		assertEquals(1, config2.processorsAfter(TestPhases.PHASE_1).size());
		assertEquals(1, config2.processorsBefore(TestPhases.PHASE_2).size());
		assertEquals(2, config2.processorsAfter(TestPhases.PHASE_2).size());
	}

}
