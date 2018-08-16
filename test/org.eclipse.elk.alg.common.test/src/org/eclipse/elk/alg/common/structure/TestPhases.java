/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * which accompanies this distribution, and is available at
 * http://www.eclipse.org/legal/epl-v10.html
 *******************************************************************************/
package org.eclipse.elk.alg.common.structure;

import org.eclipse.elk.core.util.IElkProgressMonitor;

/**
 * Phases enumeration to be used for tests.
 */
public enum TestPhases implements ILayoutPhaseFactory<TestPhases, StringBuffer> {
	PHASE_1,
	PHASE_2;

    /* (non-Javadoc)
     * @see org.eclipse.elk.core.alg.ILayoutPhaseFactory#create()
     */
    @Override
    public ILayoutPhase<TestPhases, StringBuffer> create() {
        return new ILayoutPhase<TestPhases, StringBuffer>() {
            @Override
            public void process(final StringBuffer graph, final IElkProgressMonitor progressMonitor) {
                graph.append(TestPhases.this.toString());
            }

            @Override
            public LayoutProcessorConfiguration<TestPhases, StringBuffer> getLayoutProcessorConfiguration(
                    final StringBuffer graph) {
                
                LayoutProcessorConfiguration<TestPhases, StringBuffer> config = LayoutProcessorConfiguration.create();

                // PHASE_1 -> PROCESSOR_1, PHASE_2 -> PROCESSOR_2
                config.addBefore(TestPhases.this, TestProcessors.class.getEnumConstants()[TestPhases.this.ordinal()]);
                return config;
            }
        };
    }
}
