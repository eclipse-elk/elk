/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.graph.json.test

import com.google.gson.JsonParser
import org.eclipse.elk.core.RecursiveGraphLayoutEngine
import org.eclipse.elk.core.util.BasicProgressMonitor
import org.eclipse.elk.core.util.Maybe
import org.eclipse.elk.graph.json.ElkGraphJson
import org.eclipse.elk.graph.json.JsonImporter
import org.junit.Test

import static org.junit.Assert.*
import com.google.gson.Gson

/**
 */
class SectionsTest {
    
    @Test
    def void edgeSectionTest() {
        val graph = '''
        {
          "id": "root",
          "properties": {
            "direction": "DOWN"
          },
          "children": [{"id": "n1"}, {"id": "n2"}],
          "edges": [{"id": "e1", 
            "sources": [ "n1" ], 
            "targets": [ "n2" ],
            "sections": [{
                "id": "s1",
                "startPoint": {x: 1, y: 1},
                "endPoint": {x: 2, y: 2},
                "incomingShape": "n1",
                "outgoingSections": [ 44 ]
            },{
                "id": 44,
                "startPoint": {x: 3, y: 3},
                "endPoint": {x: 4, y: 4},
                "incomingSections": [ "s1" ],
                "outgoingShape": "n2"
            }]
          }]
        }
        '''
        
        val root = ElkGraphJson.forGraph(graph).toElk
        assertTrue(root.children.size === 2)
        val n1 = root.children.findFirst[ identifier == "n1" ]
        val n2 = root.children.findFirst[ identifier == "n2" ]        
        
        assertTrue(root.containedEdges.size === 1)
        val edge = root.containedEdges.head
        assertTrue(edge.sections.size === 2)
        
        val s1 = edge.sections.findFirst[ identifier == "s1" ]
        val s44 = edge.sections.findFirst[ identifier == "44" ]
        
        assertTrue(s1.incomingShape == n1)
        assertTrue(s44.outgoingShape == n2)
        
        assertTrue(s1.outgoingSections.size === 1)
        assertTrue(s1.outgoingSections.head === s44)
        
        assertTrue(s44.incomingSections.size === 1)
        assertTrue(s44.incomingSections.head === s1)
    }
    
    @Test
    def void preserveSectionIdsTest() {
        val graph = '''
        {
          id: "root",
          properties: { 
            'elk.algorithm': 'random'
          },
          children: [
            { id: "n1", width: 10, height: 10 },
            { id: "n2", width: 10, height: 10 }
          ],
          edges: [{
            id: "e1", sources: [ "n1" ], targets: [ "n2" ],
            sections: [{
              id: "xyz",
              startPoint: { x: 0, y: 0 },
              bendPoints: [{ x: 20, y: 0 }],
              endPoint: { x: 50, y: 0 }
            }]
          }]
        }'''
        
        
        // string -> json
        val json = JsonParser.parseString(graph).asJsonObject
        
        // json -> elk
        val importer = new Maybe<JsonImporter>()
        val root = ElkGraphJson.forGraph(json).rememberImporter(importer).toElk
        val edge = root.containedEdges.head
        val es1 = edge.sections.findFirst[ identifier == "xyz" ]
        assertNotNull(es1)

        // layout
        new RecursiveGraphLayoutEngine().layout(root, new BasicProgressMonitor)
        
        importer.get.transferLayout(root)
        
        val s = new Gson().toJson(json)
        assertTrue(s.contains('''"id":"xyz"'''))
    }
    
    @Test
    def void doNotAddEmptySectionsArrayExtendedEdge() {
        val graph = '''
        {
          id: "root",
          children: [
            { id: "n1", width: 10, height: 10 },
            { id: "n2", width: 10, height: 10 }
          ],
          edges: [{ id: "e1", sources: [ "n1" ], targets: [ "n2" ] }]
        }'''
        
        val s = graph.toJsonGraphAndBackToString
        
        assertTrue(!s.contains('''"sections"'''))
    }
    
    @Test
    def void doNotAddEmptySectionsArrayPrimitiveEdge() {
        val graph = '''
        {
          id: "root",
          children: [
            { id: "n1", width: 10, height: 10 },
            { id: "n2", width: 10, height: 10 }
          ],
          edges: [{ id: "e1", source: "n1", target: "n2" }]
        }'''
        
        val s = graph.toJsonGraphAndBackToString
        
        assertTrue(!s.contains('''"sections"'''))
    }
    
    /**
     * Convert the string to an ElkNode, call JsonImporter#transferLayout 
     * and turn the ElkNode into a json string again. 
     */
    private def toJsonGraphAndBackToString(String graph) {
        // string -> json
        val json = JsonParser.parseString(graph).asJsonObject
        
        // json -> elk
        val importer = new Maybe<JsonImporter>()
        val root = ElkGraphJson.forGraph(json).rememberImporter(importer).toElk
        // transfer
        importer.get.transferLayout(root)
        // to string
        val s = new Gson().toJson(json)
        return s
    }
    
}