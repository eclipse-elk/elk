/*******************************************************************************
 * Copyright (c) 2023 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later 
 *******************************************************************************/
package org.eclipse.elk.alg.radial.intermediate.rotation;

import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.graph.ElkNode;

/**
 * The layout processor for rotation. Sets up logging and calls the angle rotation implementation.
 *
 */
public class GeneralRotator implements ILayoutProcessor<ElkNode> {
    
    @Override
    public void process(final ElkNode graph, final IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("General 'Rotator", 1);
        // elkjs-exclude-start
        if (progressMonitor.isLoggingEnabled()) {
            progressMonitor.logGraph(graph, "Before");
        }
        // elkjs-exclude-end
        IRadialRotator rotator = new AngleRotation();
        rotator.rotate(graph);
        // elkjs-exclude-start
        if (progressMonitor.isLoggingEnabled()) {
            progressMonitor.logGraph(graph, "After");
        }
        // elkjs-exclude-end
    }

}
