/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.graph.json.test

import org.eclipse.elk.alg.test.PlainJavaInitialization
import org.eclipse.elk.core.RecursiveGraphLayoutEngine
import org.eclipse.elk.core.options.CoreOptions
import org.eclipse.elk.core.util.IndividualSpacings
import org.eclipse.elk.core.util.NullElkProgressMonitor
import org.eclipse.elk.graph.json.ElkGraphJson
import org.junit.BeforeClass
import org.junit.Test

import static org.eclipse.elk.graph.util.ElkGraphUtil.*
import static org.junit.Assert.*

/**
 * Tests the proper im-/exporting of the somewhat special {@link IndividualSpacings} layout option.
 */
class IndividualSpacingsTest {
    
    static val DOUBLE_EQ = 0.001d
        
    @BeforeClass
    static def void init() {
        PlainJavaInitialization.initializePlainJavaLayout
    }
    
    @Test
    def void importIndividualSpacingsTest() {
        val graph = '''
        {
          "id": "n0",
          "children": [
            {
              "id": "outer",
              "layoutOptions": {
                  "nodeLabels.padding": "[top=0.0,left=0.0,bottom=0.0,right=0.0]"
              },
              "children": [
                {
                  "id": "i1",
                  "labels": [
                    {
                      "text": "Node 1",
                      "width": 40.0,
                      "height": 15.0
                    }
                  ],
                  "layoutOptions": {
                    "nodeLabels.placement": "[H_CENTER, V_TOP, INSIDE]"
                  },
                  "width": 60.0,
                  "height": 40.0
                },
                {
                  "id": "i2",
                  "layoutOptions": {
                    "nodeLabels.placement": "[H_CENTER, V_TOP, INSIDE]"
                  },
                  "individualSpacings": {
                    "nodeLabels.padding": "[top=10.0,left=0.0,bottom=0.0,right=0.0]"
                  },
                  "labels": [
                    {
                      "text": "Node 2",
                      "width": 40.0,
                      "height": 15.0
                    }
                  ],
                  "width": 60.0,
                  "height": 40.0
                }
              ]
            }
          ]
        }
        '''
        
        val root = ElkGraphJson.forGraph(graph).toElk
        assertTrue(root.children.size == 1) 
        val outer = root.children.findFirst[ identifier == "outer" ]
        assertTrue(outer.children.size === 2)
        val i1 = outer.children.findFirst[ identifier == "i1" ]
        val i2 = outer.children.findFirst[ identifier == "i2" ]        
        
        assertFalse(i1.hasProperty(CoreOptions.SPACING_INDIVIDUAL))
        assertTrue(i2.hasProperty(CoreOptions.SPACING_INDIVIDUAL))
        
        new RecursiveGraphLayoutEngine().layout(root, new NullElkProgressMonitor())
        val l1 = i1.labels.head
        val l2 = i2.labels.head
        
        assertEquals(0.0, l1.y, DOUBLE_EQ)
        assertEquals(10.0, l2.y, DOUBLE_EQ)
    }
    
    @Test
    def void exportIndividualSpacings() {
        val graph = createGraph

        graph.setProperty(CoreOptions.SPACING_NODE_NODE, 10.0)

        val individualSpacings = new IndividualSpacings()
        individualSpacings.setProperty(CoreOptions.SPACING_NODE_NODE, 20.0)
        graph.setProperty(CoreOptions.SPACING_INDIVIDUAL, individualSpacings)

        val json = ElkGraphJson.forGraph(graph)
                               .omitUnknownLayoutOptions(true)
                               .toJson

        assertTrue(json.contains("individualSpacings"))
        assertTrue(json.contains("nodeNode\":\"10"))
        assertTrue(json.contains("nodeNode\":\"20"))
    }
    
    @Test
    def exportNoIndividualSpacings() {
        val graph = createGraph
        val json = ElkGraphJson.forGraph(graph)
                              .omitUnknownLayoutOptions(true)
                              .toJson
        assertFalse(json.contains("individualSpacings"))
        assertFalse(json.contains("IndividualSpacings"))
    }
    
}