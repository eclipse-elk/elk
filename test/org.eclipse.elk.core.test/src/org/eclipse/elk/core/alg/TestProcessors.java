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

import org.eclipse.elk.core.util.IElkProgressMonitor;

public enum TestProcessors implements ILayoutProcessorFactory<StringBuffer> {
	PROCESSOR_1,
	PROCESSOR_2,
	PROCESSOR_3;
	
	public ILayoutProcessor<StringBuffer> create() {
		return new ILayoutProcessor<StringBuffer>() {
			@Override
			public void process(StringBuffer graph, IElkProgressMonitor progressMonitor) {
			    graph.append(TestProcessors.this.toString());
			}
		};
	}
}
