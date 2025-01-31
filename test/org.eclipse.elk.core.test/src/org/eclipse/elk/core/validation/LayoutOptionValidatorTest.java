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

import org.eclipse.elk.core.data.LayoutMetaDataService;
import org.eclipse.elk.core.data.LayoutOptionData;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.util.ElkUtil;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.util.ElkGraphUtil;
import org.junit.Test;

/**
 * Tests for {@link LayoutOptionValidator}.
 */
public class LayoutOptionValidatorTest {
    
    @Test
    public void testInvalidType() {
        ElkNode graph = ElkGraphUtil.createGraph();
        LayoutOptionData aspectRatio = LayoutMetaDataService.getInstance().getOptionData(CoreOptions.ASPECT_RATIO.getId());
        graph.setProperty(aspectRatio, "foo");
        
        LayoutOptionValidator validator = new LayoutOptionValidator();
        ElkUtil.applyVisitors(graph, validator);
        
        assertEquals(
                "[ERROR: The assigned value foo of the option 'Aspect Ratio' does not match the type Double. (Root Node)]",
                validator.getIssues().toString());
    }
    
    @Test
    public void testExclusiveLowerBound() {
        ElkNode graph = ElkGraphUtil.createGraph();
        graph.setProperty(CoreOptions.ASPECT_RATIO, 0.0);
        
        LayoutOptionValidator validator = new LayoutOptionValidator();
        ElkUtil.applyVisitors(graph, validator);
        
        assertEquals(
                "[ERROR: The assigned value 0.0 of the option 'org.eclipse.elk.aspectRatio' is less than the lower bound 0.0 (exclusive). (Root Node)]",
                validator.getIssues().toString());
    }

}
