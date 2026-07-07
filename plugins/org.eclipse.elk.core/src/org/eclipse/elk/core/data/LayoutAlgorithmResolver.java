/*******************************************************************************
 * Copyright (c) 2018, 2020 TypeFox GmbH and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.data;

import org.eclipse.elk.core.UnsupportedConfigurationException;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.util.ElkUtil;
import org.eclipse.elk.core.util.IGraphElementVisitor;
import org.eclipse.elk.graph.ElkGraphElement;
import org.eclipse.elk.graph.ElkNode;

/**
 * Resolves layout algorithms configured through the {@link CoreOptions#ALGORITHM} option and assigns
 * the resulting meta data to the {@link CoreOptions#RESOLVED_ALGORITHM} option.
 */
public class LayoutAlgorithmResolver implements IGraphElementVisitor {

    @Override
    public void visit(final ElkGraphElement element) {
        if (element instanceof ElkNode && !element.getProperty(CoreOptions.NO_LAYOUT)) {
            resolveAlgorithm((ElkNode) element);
        }
    }
    
    /**
     * Resolve the layout algorithm to apply to the content of the given node. If no algorithm is configured for the
     * given node but should be, we try to resolve the {@link #getDefaultLayoutAlgorithmID() default algorithm}. If the
     * algorithm configured for the node does not exist, an exception is thrown instead.
     * 
     * @param node
     *            the node to resolve the layout algorithm for.
     * @throws UnsupportedConfigurationException
     *             if a non-existent layout algorithm was configured for the node.
     */
    protected void resolveAlgorithm(final ElkNode node) {
        String algorithmId = node.getProperty(CoreOptions.ALGORITHM);
        
        // Stage 1: Try to resolve the intended algorithm
        if (resolveAndSetAlgorithm(algorithmId, node)) {
            return;
        }
        
        // Stage 2: If we must resolve a layout algorithm, try to fall back on the default if none was specified
        if (mustResolve(node)) {
             if (algorithmId == null || algorithmId.trim().isEmpty()) {
                 // No algorithm was specified, load the default one
                 String defaultAlgorithmId = getDefaultLayoutAlgorithmID();
                 
                 if (!resolveAndSetAlgorithm(getDefaultLayoutAlgorithmID(), node)) {
                     StringBuilder message = new StringBuilder("Unable to load default layout algorithm ")
                             .append(defaultAlgorithmId)
                             .append(" for unconfigured node ");
                     ElkUtil.printElementPath(node, message);
                     throw new UnsupportedConfigurationException(message.toString());
                 }
             
             } else {
                 // An algorithm was specified, but not found. Fail!
                 StringBuilder message = new StringBuilder("Layout algorithm '")
                         .append(algorithmId)
                         .append("' not found for ");
                 ElkUtil.printElementPath(node, message);
                 throw new UnsupportedConfigurationException(message.toString());
             }
        }
    }
    
    /**
     * Tries to resolve the layout algorithm with the given ID. If one is found, the node's
     * {@link CoreOptions#RESOLVED_ALGORITHM} is set accordingly.
     * 
     * @param algorithmId
     *            the algorithm to resolve.
     * @param node
     *            the node to resolve it for.
     * @return {@code true} if the algorithm was resolved and stored.
     */
    protected boolean resolveAndSetAlgorithm(final String algorithmId, final ElkNode node) {
        LayoutAlgorithmData algorithmData = LayoutMetaDataService.getInstance().getAlgorithmDataBySuffix(algorithmId);
        
        if (algorithmData != null) {
            node.setProperty(CoreOptions.RESOLVED_ALGORITHM, algorithmData);
            return true;
        } else {
            return false;
        }
    }

    /**
     * Determine whether it is mandatory to resolve the layout algorithm for the given node.
     */
    protected boolean mustResolve(final ElkNode node) {
        // Hierarchical nodes must have their layout algorithm resolved, but so do nodes that have inside self loops
        // self loop and nothing else.
        return !node.hasProperty(CoreOptions.RESOLVED_ALGORITHM)
                && (!node.getChildren().isEmpty() || node.getProperty(CoreOptions.INSIDE_SELF_LOOPS_ACTIVATE));
    }
    
    /**
     * Returns the ID of the layout algorithm to be used if no algorithm is explicitly set in the graph.
     * 
     * <p>The default implementation returns {@code "org.eclipse.elk.layered"}.</p>
     * 
     * @return ID of the default layout algorithm.
     */
    public String getDefaultLayoutAlgorithmID() {
        return "org.eclipse.elk.layered";
    }
    
}
