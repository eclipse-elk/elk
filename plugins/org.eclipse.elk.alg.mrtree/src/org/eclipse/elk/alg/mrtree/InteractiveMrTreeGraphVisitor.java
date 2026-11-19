/*******************************************************************************
 * Copyright (c) 2022 Kiel University and others.
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
package org.eclipse.elk.alg.mrtree;

import org.eclipse.elk.alg.mrtree.options.MrTreeOptions;
import org.eclipse.elk.alg.mrtree.options.OrderWeighting;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.util.IGraphElementVisitor;
import org.eclipse.elk.graph.ElkGraphElement;
import org.eclipse.elk.graph.ElkNode;

/**
 * Graph visitor which visits only the root node and recursively steps through the graph
 * to set interactive options for the mrtree algorithm.
 * The mrtree algorithm needs {@code Interactive} to be set to true in all root nodes.
 */
public class InteractiveMrTreeGraphVisitor implements IGraphElementVisitor {

    /* (non-Javadoc)
     * @see org.eclipse.elk.core.util.IGraphElementVisitor#visit(org.eclipse.elk.graph.ElkGraphElement)
     */
    @Override
    public void visit(final ElkGraphElement element) {
        // Only apply to root of the graph
        if (element instanceof ElkNode) {
            ElkNode root = (ElkNode) element;
            setInteractiveOptions(root);
        }
    }
    
    /**
     * Sets the required options for the interactive layout run with the {@code rectpacking} algorithm.
     */
    public void setInteractiveOptions(final ElkNode root) {
        String algorithm = root.getProperty(CoreOptions.ALGORITHM);
        if (algorithm != null && MrTreeOptions.ALGORITHM_ID.endsWith(algorithm) && !root.getChildren().isEmpty()) {
            root.setProperty(CoreOptions.INTERACTIVE, true);
            root.setProperty(MrTreeOptions.WEIGHTING, OrderWeighting.CONSTRAINT);
        }
        return;
    }

}
