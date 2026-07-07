/*******************************************************************************
 * Copyright (c) 2018 TypeFox GmbH and others.
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
package org.eclipse.elk.core.validation;

import java.lang.reflect.InvocationTargetException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.eclipse.elk.core.data.LayoutAlgorithmData;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.util.ElkUtil;
import org.eclipse.elk.graph.ElkConnectableShape;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkEdgeSection;
import org.eclipse.elk.graph.ElkGraphElement;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.util.ElkGraphUtil;

/**
 * A graph validator that checks basic structural properties of a graph and applies algorithm-specific
 * validation.
 */
public class GraphValidator implements IValidatingGraphElementVisitor {
    
    /** The issues found by this validator. */
    private final List<GraphIssue> issues = new ArrayList<GraphIssue>();
    /** A cache of algorithm-specific validator instances. */
    private final Map<LayoutAlgorithmData, IValidatingGraphElementVisitor> algorithmSpecificValidators
            = new HashMap<>();

    /**
     * Check the structural properties of the given graph element.
     */
    @Override
    public void visit(final ElkGraphElement element) {
        // Execute the generic checks
        if (element instanceof ElkEdge) {
            checkEdge((ElkEdge) element);
        }
        
        // Execute the algorithm-specific checks
        ElkNode parent = ElkGraphUtil.containingGraph(element);
        if (parent != null) {
            runAlgorithmSpecificChecks(element, parent);
        }
        if (element instanceof ElkNode && (parent == null || parent.getProperty(CoreOptions.RESOLVED_ALGORITHM)
                != element.getProperty(CoreOptions.RESOLVED_ALGORITHM))) {
            runAlgorithmSpecificChecks(element, (ElkNode) element);
        }
    }
    
    /**
     * Apply generic checks to the given edge.
     */
    protected void checkEdge(final ElkEdge edge) {
        if (!edge.isConnected()) {
            issues.add(new GraphIssue(edge, "Edge is not connected.", GraphIssue.Severity.ERROR));
        } else {
            ElkNode bestContainer = ElkGraphUtil.findBestEdgeContainment(edge);
            if (bestContainer != null && bestContainer != edge.getContainingNode()) {
                StringBuilder message = new StringBuilder("Edge should be contained in ");
                ElkUtil.printElementPath(bestContainer, message);
                issues.add(new GraphIssue(edge, message.toString(), GraphIssue.Severity.WARNING));
            }
        }
        
        for (ElkEdgeSection edgeSection : edge.getSections()) {
            ElkConnectableShape incomingShape = edgeSection.getIncomingShape();
            if (incomingShape != null) {
                if (!edge.getSources().contains(incomingShape)) {
                    String message = incomingShape.eClass().getName()
                            + " declared as incoming shape is not a source of this edge.";
                    issues.add(new GraphIssue(edge, message, GraphIssue.Severity.ERROR));
                }
                if (!edgeSection.getIncomingSections().isEmpty()) {
                    String message = "An edge section cannot be connected to an "
                            + incomingShape.eClass().getName() + " and other sections at the same time.";
                    issues.add(new GraphIssue(edge, message, GraphIssue.Severity.ERROR));
                }
            }
            ElkConnectableShape outgoingShape = edgeSection.getOutgoingShape();
            if (outgoingShape != null) {
                if (!edge.getTargets().contains(outgoingShape)) {
                    String message = outgoingShape.eClass().getName()
                            + " declared as outgoing shape is not a target of this edge.";
                    issues.add(new GraphIssue(edge, message, GraphIssue.Severity.ERROR));
                }
                if (!edgeSection.getOutgoingSections().isEmpty()) {
                    String message = "An edge section cannot be connected to an "
                            + outgoingShape.eClass().getName() + " and other sections at the same time.";
                    issues.add(new GraphIssue(edge, message, GraphIssue.Severity.ERROR));
                }
            }
        }
    }
    
    /**
     * Apply algorithm-specific checks to the given graph element.
     */
    protected void runAlgorithmSpecificChecks(final ElkGraphElement element, final ElkNode parent) {
        LayoutAlgorithmData algoData = parent.getProperty(CoreOptions.RESOLVED_ALGORITHM);
        if (algoData != null) {
            IValidatingGraphElementVisitor validator = getValidator(algoData);
            if (validator != null) {
                validator.visit(element);
            }
        }
    }
    
    /**
     * Create an instance of the algorithm-specific validator for the given layout algorithm.
     * The instance is cached for performance optimization.
     */
    protected IValidatingGraphElementVisitor getValidator(final LayoutAlgorithmData algoData) {
        IValidatingGraphElementVisitor validator = algorithmSpecificValidators.get(algoData);
        if (validator != null) {
            return validator;
        }
        Class<? extends IValidatingGraphElementVisitor> validatorClass = algoData.getValidatorClass();
        if (validatorClass == null) {
            return null;
        }
        try {
            validator = validatorClass.getConstructor().newInstance();
            algorithmSpecificValidators.put(algoData, validator);
            return validator;
        } catch (InstantiationException | IllegalAccessException | InvocationTargetException
                | NoSuchMethodException ex) {
            throw new RuntimeException("Failed to instantiate validator for " + algoData.getId(), ex);
        }
    }
    
    @Override
    public List<GraphIssue> getIssues() {
        if (algorithmSpecificValidators.isEmpty()) {
            return issues;
        } else {
            List<GraphIssue> result = new ArrayList<>();
            result.addAll(issues);
            for (IValidatingGraphElementVisitor validator : algorithmSpecificValidators.values()) {
                result.addAll(validator.getIssues());
            }
            return result;
        }
    }

}
