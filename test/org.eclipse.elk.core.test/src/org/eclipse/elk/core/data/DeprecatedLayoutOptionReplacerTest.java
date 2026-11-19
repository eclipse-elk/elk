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

import java.util.EnumSet;

import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.options.PortLabelPlacement;
import org.eclipse.elk.core.options.SizeOptions;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.junit.Before;
import org.junit.Test;

/**
 * Tests {@link DeprecatedLayoutOptionReplacer}.
 */
public class DeprecatedLayoutOptionReplacerTest {

    private final DeprecatedLayoutOptionReplacer replacer = new DeprecatedLayoutOptionReplacer();

    @Before
    public void initMetaDataService() {
        LayoutMetaDataService.getInstance().registerLayoutMetaDataProviders(new CoreOptions());
    }

    @SuppressWarnings("deprecation")
    @Test
    public void testReplacesNextToPortIfPossible() {
        ElkNode node = ElkGraphUtil.createGraph();
        node.setProperty(CoreOptions.PORT_LABELS_NEXT_TO_PORT_IF_POSSIBLE, true);

        replacer.visit(node);

        assertFalse(node.hasProperty(CoreOptions.PORT_LABELS_NEXT_TO_PORT_IF_POSSIBLE));
        assertTrue(node.hasProperty(CoreOptions.PORT_LABELS_PLACEMENT));
        assertTrue(node.getProperty(CoreOptions.PORT_LABELS_PLACEMENT)
                .contains(PortLabelPlacement.NEXT_TO_PORT_IF_POSSIBLE));
    }

    @SuppressWarnings("deprecation")
    @Test
    public void testMovesSpaceEfficientPortLabels() {
        ElkNode node = ElkGraphUtil.createGraph();
        node.setProperty(CoreOptions.NODE_SIZE_OPTIONS, EnumSet.of(SizeOptions.SPACE_EFFICIENT_PORT_LABELS));

        replacer.visit(node);

        assertFalse(node.getProperty(CoreOptions.NODE_SIZE_OPTIONS).contains(SizeOptions.SPACE_EFFICIENT_PORT_LABELS));
        assertTrue(node.hasProperty(CoreOptions.PORT_LABELS_PLACEMENT));
        assertTrue(node.getProperty(CoreOptions.PORT_LABELS_PLACEMENT)
                .contains(PortLabelPlacement.SPACE_EFFICIENT));        
    }
}
