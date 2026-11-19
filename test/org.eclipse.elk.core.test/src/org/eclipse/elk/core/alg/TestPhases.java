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
package org.eclipse.elk.core.alg;

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
