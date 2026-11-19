/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
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
package org.eclipse.elk.core.data;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import org.eclipse.elk.core.data.LayoutDataContentAssist.Proposal;
import org.eclipse.elk.core.options.BoxLayouterOptions;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.options.FixedLayouterOptions;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkGraphElement;
import org.eclipse.elk.graph.ElkGraphFactory;
import org.eclipse.elk.graph.ElkGraphPackage;
import org.eclipse.elk.graph.ElkLabel;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.ElkPort;
import org.eclipse.elk.graph.properties.IProperty;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.eclipse.emf.ecore.EClass;
import org.junit.Before;
import org.junit.Test;
import org.junit.experimental.runners.Enclosed;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameter;
import org.junit.runners.Parameterized.Parameters;

/**
 * Tests {@link LayoutDataContentAssist} methods.
 */
@RunWith(Enclosed.class)
public class LayoutDataContentAssistTest {

    @Before
    public void initMetaDataService() {
        LayoutMetaDataService.getInstance().registerLayoutMetaDataProviders(new CoreOptions());
    }

    @RunWith(Parameterized.class)
    public static class LayoutOptionsPrefixes {

        @Parameters(name = "Prefix: {2}")
        public static Collection<Object[]> data() {

            // @formatter:off
            return Arrays.asList(new Object[][] { 
                // Prefixes of the id
                { ElkNode.class, CoreOptions.PORT_CONSTRAINTS, "org.eclipse.elk.portConstraints" },
                { ElkNode.class, CoreOptions.PORT_CONSTRAINTS, "eclipse.elk.portConstraints" },
                { ElkNode.class, CoreOptions.PORT_CONSTRAINTS, "elk.portConstraints" },
                { ElkNode.class, CoreOptions.PORT_CONSTRAINTS, "portConstraints" },
                { ElkNode.class, CoreOptions.PORT_CONSTRAINTS, "portCo" },
                
                // Some part of the name
                { ElkNode.class, CoreOptions.PORT_CONSTRAINTS, "rt Constr" },
                
                { ElkEdge.class, CoreOptions.PRIORITY, "prio" },
                
                { ElkPort.class, CoreOptions.PORT_SIDE, "port.sid" },
                
                { ElkLabel.class, CoreOptions.NODE_LABELS_PLACEMENT, "nodeLabels.pla" },
                
                // Anything that is assinged to 'nodes' should be returned for an empty prefix 
                { ElkNode.class, CoreOptions.PADDING, "" },
                { ElkNode.class, CoreOptions.ALIGNMENT, "" },
            });
            // @formatter:on
        }

        @Parameter(0)
        public Class<?> graphElementClass;
        @Parameter(1)
        public IProperty<?> expectedOption;
        @Parameter(2)
        public String prefix;

        @Test
        public void testLayoutOption() {
            EClass clazz = (EClass) ElkGraphPackage.eINSTANCE.getEClassifier(graphElementClass.getSimpleName());
            ElkGraphElement graphElement = (ElkGraphElement) ElkGraphFactory.eINSTANCE.create(clazz);
            List<Proposal<LayoutOptionData>> ps =
                    LayoutDataContentAssist.getLayoutOptionProposals(graphElement, prefix);
            boolean contains = ps.stream().anyMatch(p -> p.data.equals(expectedOption));
            assertTrue(contains);
        }
        
    }

    public static class Other {
        
        @Test
        public void testParentNodes() {
            ElkNode root = ElkGraphUtil.createGraph();
            ElkNode child = ElkGraphUtil.createNode(root);

            List<Proposal<LayoutOptionData>> parentPs =
                    LayoutDataContentAssist.getLayoutOptionProposals(root, CoreOptions.ALGORITHM.getId());
            boolean parentContains = parentPs.stream().anyMatch(p -> p.data.equals(CoreOptions.ALGORITHM));
            assertTrue(parentContains);

            List<Proposal<LayoutOptionData>> childPs =
                    LayoutDataContentAssist.getLayoutOptionProposals(child, CoreOptions.ALGORITHM.getId());
            boolean childContains = childPs.stream().anyMatch(p -> p.data.equals(CoreOptions.ALGORITHM));
            assertFalse(childContains);
        }

        @Test
        public void testAlgorithmSpecificKnown() {
            ElkNode root = ElkGraphUtil.createGraph();
            ElkGraphUtil.createNode(root);
            root.setProperty(CoreOptions.ALGORITHM, BoxLayouterOptions.ALGORITHM_ID);

            List<Proposal<LayoutOptionData>> parentPs =
                    LayoutDataContentAssist.getLayoutOptionProposals(root, CoreOptions.EXPAND_NODES.getId());
            boolean parentContains = parentPs.stream().anyMatch(p -> p.data.equals(CoreOptions.EXPAND_NODES));
            assertTrue(parentContains);
        }

        @Test
        public void testAlgorithmSpecificUnknown() {
            ElkNode root = ElkGraphUtil.createGraph();
            ElkGraphUtil.createNode(root);
            root.setProperty(CoreOptions.ALGORITHM, FixedLayouterOptions.ALGORITHM_ID);

            List<Proposal<LayoutOptionData>> parentPs =
                    LayoutDataContentAssist.getLayoutOptionProposals(root, CoreOptions.EXPAND_NODES.getId());
            boolean parentContains = parentPs.stream().anyMatch(p -> p.data.equals(CoreOptions.EXPAND_NODES));
            assertFalse(parentContains);
        }
        
        @Test
        public void testOptionGroupIsAdded() {
            ElkNode root = ElkGraphUtil.createGraph();
            ElkGraphUtil.createNode(root);
            List<Proposal<LayoutOptionData>> ps = LayoutDataContentAssist.getLayoutOptionProposals(root, "nodeN");
            // The option must be included _and_ the proposal must contain the group
            boolean groupIsAdded = ps.stream().anyMatch(
                    p -> p.data.equals(CoreOptions.SPACING_NODE_NODE) && p.proposal.contains(p.data.getGroup()));
            assertTrue(groupIsAdded);
        }
        
    }
}
