/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.graph.json.test

import org.eclipse.elk.graph.json.ElkGraphJson
import org.junit.Test

import static org.junit.Assert.*

/**
 */
class EdgesTest {
    
    @Test
    def void extendedEdgeTest() {
        val graph = '''
        {
          "id": "root",
          "children": [{"id": "n1"}, {"id": 3}],
          "edges": [{"id": "e1", 
            "sources": [ "n1" ], 
            "targets": [ 3 ]
          }]
        }
        '''
        
        val root = ElkGraphJson.forGraph(graph).toElk
        assertTrue(root.children.size === 2)
        assertTrue(root.containedEdges.size === 1)

        val n1 = root.children.findFirst[c | c.identifier == "n1" ]
        val n3 = root.children.findFirst[c | c.identifier == "3" ]
        
        val edge = root.containedEdges.head
        assertTrue(edge.sources.size === 1)
        assertTrue(edge.sources.head === n1)
        assertTrue(edge.targets.size === 1)
        assertTrue(edge.targets.head === n3)        
    }
    
    @Test 
    def void primitiveEdgeTest() {
        val graph = '''
        {
          "id": "root",
          "children": [{"id": "n1"}, {"id": 3}],
          "edges": [{"id": "e1", 
            "source": "n1", 
            "target": 3
          }]
        }
        '''
        
        val root = ElkGraphJson.forGraph(graph).toElk
        assertTrue(root.children.size === 2)
        assertTrue(root.containedEdges.size === 1)

        val n1 = root.children.findFirst[c | c.identifier == "n1" ]
        val n3 = root.children.findFirst[c | c.identifier == "3" ]
        
        val edge = root.containedEdges.head
        assertTrue(edge.sources.size === 1)
        assertTrue(edge.sources.head === n1)
        assertTrue(edge.targets.size === 1)
        assertTrue(edge.targets.head === n3)        
    }

    @Test
    def void edgeContainmentTest() {
        val graph = '''
        {
          "id": "root",
          "children": [
            {"id": "p"},
            {"id": "q",
              "children": [{"id": "r"}],
              "edges": [{"id": "e", "source": "p", "target": "r" }]
            }]
        }
        '''
        
        val root = ElkGraphJson.forGraph(graph).toElk
        assertTrue(root.containedEdges.size === 1)
    }
}
