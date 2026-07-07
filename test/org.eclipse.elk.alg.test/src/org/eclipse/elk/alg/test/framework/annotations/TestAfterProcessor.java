/*******************************************************************************
 * Copyright (c) 2018, 2019 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.test.framework.annotations;

import java.lang.annotation.ElementType;
import java.lang.annotation.Repeatable;
import java.lang.annotation.Retention;
import java.lang.annotation.RetentionPolicy;
import java.lang.annotation.Target;

import org.eclipse.elk.core.alg.ILayoutProcessor;

/**
 * Flags a method as a whitebox test. The test will be run whenever the given processor just finished executing. This
 * annotation can be used more then once per class or method. A method annotated with this annotation must expect
 * exactly one parameter of type {@code Object}.
 */
@Retention(RetentionPolicy.RUNTIME)
@Target(ElementType.METHOD)
@Repeatable(TestAfterProcessors.class)
public @interface TestAfterProcessor {
    
    /** The processor the graph should be executed after. */
    Class<? extends ILayoutProcessor<?>> value();
    
}
