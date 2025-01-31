/*******************************************************************************
 * Copyright (c) 2018 TypeFox GmbH and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.core.validation;

import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import org.eclipse.elk.core.AbstractLayoutProvider;
import org.eclipse.elk.core.data.LayoutAlgorithmData;
import org.eclipse.elk.core.data.LayoutAlgorithmResolver;
import org.eclipse.elk.core.data.LayoutMetaDataService;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.util.AlgorithmFactory;
import org.eclipse.elk.core.util.ElkUtil;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.graph.ElkGraphElement;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.junit.Test;

/**
 * Tests for {@link GraphValidator}.
 */
public class GraphValidatorTest {
    
    @Test
    public void testUnconnectedEdge() {
        ElkNode graph = ElkGraphUtil.createGraph();
        ElkGraphUtil.createEdge(graph);
        
        GraphValidator validator = new GraphValidator();
        ElkUtil.applyVisitors(graph, validator);
        
        assertEquals("[ERROR: Edge is not connected. (Root Node > Edge)]", validator.getIssues().toString());
    }
    
    @Test
    public void testCustomValidation() {
        LayoutMetaDataService.getInstance().registerLayoutMetaDataProviders(registry -> {
            registry.register(
                new LayoutAlgorithmData.Builder()
                    .id("foo.algorithm")
                    .providerFactory(new AlgorithmFactory(FooAlgorithm.class))
                    .validatorClass(FooValidator.class)
                    .create()
            );
        });
        
        ElkNode graph = ElkGraphUtil.createGraph();
        ElkGraphUtil.createNode(graph).setIdentifier("foo");
        graph.setProperty(CoreOptions.ALGORITHM, "foo.algorithm");
        
        LayoutAlgorithmResolver algorithmResolver = new LayoutAlgorithmResolver();
        GraphValidator validator = new GraphValidator();
        ElkUtil.applyVisitors(graph, algorithmResolver, validator);
        
        assertEquals("[ERROR: FOO! (Root Node > Node foo)]", validator.getIssues().toString());
    }
    
    public static class FooAlgorithm extends AbstractLayoutProvider {
        @Override
        public void layout(ElkNode layoutGraph, IElkProgressMonitor progressMonitor) {
            layoutGraph.setDimensions(100, 100);
        }
    }
    
    public static class FooValidator implements IValidatingGraphElementVisitor {
        private final List<GraphIssue> issues = new ArrayList<GraphIssue>();
        
        @Override
        public void visit(ElkGraphElement element) {
            if ("foo".equals(element.getIdentifier())) {
                issues.add(new GraphIssue(element, "FOO!", GraphIssue.Severity.ERROR));
            }
        }

        @Override
        public Collection<GraphIssue> getIssues() {
            return issues;
        }
    }

}
